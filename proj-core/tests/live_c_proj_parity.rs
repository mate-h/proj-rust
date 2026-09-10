#![cfg(feature = "c-proj-compat")]

#[path = "common/c_proj_ffi.rs"]
mod c_proj_ffi;

use c_proj_ffi::CProjTransform;
use proj::Proj;
use proj_core::{AreaOfInterest, Coord, SelectionOptions, Transform};
use serde::Deserialize;

#[derive(Deserialize)]
struct ReferencePoint {
    from_epsg: u32,
    to_epsg: u32,
    input_x: f64,
    input_y: f64,
    input_z: Option<f64>,
    expected_x: f64,
    expected_y: f64,
    expected_z: Option<f64>,
    tolerance: f64,
    tolerance_z: Option<f64>,
    description: String,
}

impl ReferencePoint {
    fn is_3d(&self) -> bool {
        self.input_z.is_some() && self.expected_z.is_some()
    }
}

fn load_corpus() -> Vec<ReferencePoint> {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../testdata/reference_values.json"
    );
    let data =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));
    serde_json::from_str(&data).unwrap_or_else(|e| panic!("failed to parse {path}: {e}"))
}

fn geographic_aoi(x: f64, y: f64) -> Option<(f64, f64)> {
    ((-180.0..=180.0).contains(&x) && (-90.0..=90.0).contains(&y)).then_some((x, y))
}

/// Live C PROJ result: XY, plus Z when the corpus record is 3D.
fn live_c_proj(point: &ReferencePoint) -> Result<(f64, f64, Option<f64>), String> {
    let from = format!("EPSG:{}", point.from_epsg);
    let to = format!("EPSG:{}", point.to_epsg);
    if point.is_3d() {
        let z = point.input_z.expect("3D record has input_z");
        let aoi = geographic_aoi(point.input_x, point.input_y);
        CProjTransform::new_promoted_3d(point.from_epsg, point.to_epsg, aoi)
            .and_then(|proj| proj.convert_3d((point.input_x, point.input_y, z)))
            .map(|(x, y, z)| (x, y, Some(z)))
            .map_err(|e| format!("C PROJ convert_3d failed for {from}->{to}: {e}"))
    } else {
        let proj = Proj::new_known_crs(&from, &to, None)
            .map_err(|e| format!("failed to create C PROJ transform {from}->{to}: {e}"))?;
        proj.convert((point.input_x, point.input_y))
            .map(|(x, y)| (x, y, None))
            .map_err(|e| format!("C PROJ convert failed for {from}->{to}: {e}"))
    }
}

fn rust_transform(point: &ReferencePoint) -> proj_core::Result<Transform> {
    if !point.is_3d() {
        return Transform::from_epsg(point.from_epsg, point.to_epsg);
    }

    let mut options = SelectionOptions::default();
    if let Some((x, y)) = geographic_aoi(point.input_x, point.input_y) {
        options.area_of_interest = Some(AreaOfInterest::source_crs_point(Coord::new(x, y)));
    }
    Transform::with_selection_options(
        &format!("EPSG:{}", point.from_epsg),
        &format!("EPSG:{}", point.to_epsg),
        options,
    )
}

fn rust_convert(
    transform: &Transform,
    point: &ReferencePoint,
) -> Result<(f64, f64, Option<f64>), String> {
    if point.is_3d() {
        let z = point.input_z.expect("3D record has input_z");
        transform
            .convert_3d((point.input_x, point.input_y, z))
            .map(|(x, y, z)| (x, y, Some(z)))
            .map_err(|e| format!("{}: proj-core convert_3d failed: {e}", point.description))
    } else {
        transform
            .convert((point.input_x, point.input_y))
            .map(|(x, y)| (x, y, None))
            .map_err(|e| format!("{}: proj-core convert failed: {e}", point.description))
    }
}

fn assert_within_tolerance(
    point: &ReferencePoint,
    expected: (f64, f64, Option<f64>),
    actual: (f64, f64, Option<f64>),
) -> Option<String> {
    let dx = (actual.0 - expected.0).abs();
    let dy = (actual.1 - expected.1).abs();
    let dz = match (expected.2, actual.2) {
        (Some(e), Some(a)) => Some((a - e).abs()),
        _ => None,
    };
    let tol_z = point.tolerance_z.unwrap_or(point.tolerance);
    if dx <= point.tolerance && dy <= point.tolerance && dz.is_none_or(|d| d <= tol_z) {
        return None;
    }

    Some(format!(
        "{}: expected {:?}, got {:?}, delta ({:e}, {:e}, {:?}), tol {:e}, tol_z {:e}",
        point.description,
        (expected.0, expected.1, expected.2),
        (actual.0, actual.1, actual.2),
        dx,
        dy,
        dz,
        point.tolerance,
        tol_z
    ))
}

#[test]
fn reference_corpus_stays_in_sync_with_live_c_proj() {
    let corpus = load_corpus();
    assert!(!corpus.is_empty(), "corpus is empty");

    let mut failures = Vec::new();

    for point in &corpus {
        let actual = match live_c_proj(point) {
            Ok(actual) => actual,
            Err(error) => {
                failures.push(format!("{}: {error}", point.description));
                continue;
            }
        };
        if let Some(failure) = assert_within_tolerance(
            point,
            (point.expected_x, point.expected_y, point.expected_z),
            actual,
        ) {
            failures.push(failure);
        }
    }

    if !failures.is_empty() {
        panic!(
            "{} of {} corpus points drifted from live C PROJ:\n{}",
            failures.len(),
            corpus.len(),
            failures.join("\n")
        );
    }
}

#[test]
fn proj_core_matches_live_c_proj_for_supported_corpus_cases() {
    let corpus = load_corpus();
    assert!(!corpus.is_empty(), "corpus is empty");

    let mut compared = 0;
    let mut skipped = 0;
    let mut failures = Vec::new();

    for point in &corpus {
        let transform = match rust_transform(point) {
            Ok(transform) => transform,
            Err(_) => {
                skipped += 1;
                continue;
            }
        };

        let expected = match live_c_proj(point) {
            Ok(expected) => expected,
            Err(error) => {
                failures.push(format!("{}: {error}", point.description));
                continue;
            }
        };
        let actual = match rust_convert(&transform, point) {
            Ok(actual) => actual,
            Err(error) => {
                failures.push(error);
                continue;
            }
        };

        if let Some(failure) = assert_within_tolerance(point, expected, actual) {
            failures.push(failure);
        } else {
            compared += 1;
        }
    }

    assert!(
        compared >= 100,
        "expected broad live coverage, only compared {compared} points ({skipped} skipped)"
    );

    if !failures.is_empty() {
        panic!(
            "{} of {} live C PROJ comparisons failed ({} skipped):\n{}",
            failures.len(),
            corpus.len(),
            skipped,
            failures.join("\n")
        );
    }
}
