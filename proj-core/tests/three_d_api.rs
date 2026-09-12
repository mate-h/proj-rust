use proj_core::{
    CompoundCrsDef, Coord3D, CoordinateOperationId, CrsDef, Datum, LinearUnit, OperationDomain,
    Transform, VerticalCrsDef,
};

/// WGS 84 geodetic (lon/lat/h) → ECEF checkpoint for NYC from C PROJ and
/// GeographicLib cartography on the WGS 84 ellipsoid.
const NYC_LON_LAT_H: (f64, f64, f64) = (-74.006, 40.7128, 10.0);
const NYC_ECEF_XYZ: (f64, f64, f64) = (
    1_334_000.544_686_07,
    -4_654_052.129_206_82,
    4_138_306.761_372_84,
);

#[test]
fn wgs84_geographic_to_ecef_checkpoint() {
    let t = Transform::new("EPSG:4326", "EPSG:4978").unwrap();
    let (x, y, z) = t.convert_3d(NYC_LON_LAT_H).unwrap();

    assert!((x - NYC_ECEF_XYZ.0).abs() < 1e-4, "x = {x}");
    assert!((y - NYC_ECEF_XYZ.1).abs() < 1e-4, "y = {y}");
    assert!((z - NYC_ECEF_XYZ.2).abs() < 1e-4, "z = {z}");

    let inv = t.inverse().unwrap();
    let (lon, lat, h) = inv.convert_3d((x, y, z)).unwrap();
    assert!((lon - NYC_LON_LAT_H.0).abs() < 1e-10);
    assert!((lat - NYC_LON_LAT_H.1).abs() < 1e-10);
    assert!((h - NYC_LON_LAT_H.2).abs() < 1e-6);
}

#[test]
fn cross_datum_ellipsoidal_3d_to_ecef() {
    // ETRS89 3D to WGS 84 ECEF: height is consumed by cart framing.
    let compound = Transform::new("EPSG:4937", "EPSG:4978").unwrap();
    assert!(
        compound
            .vertical_diagnostics()
            .operation_name
            .as_deref()
            .is_some_and(|name| name.contains("embedded in geocentric")),
        "expected geocentric vertical embedding, got {:?}",
        compound.vertical_diagnostics()
    );

    // Both-3D endpoints apply Helmert in 3D, so this matches the ECEF path.
    let from_compound = compound.convert_3d(NYC_LON_LAT_H).unwrap();
    let via_etrs89_ecef = {
        let to_etrs89 = Transform::new("EPSG:4937", "EPSG:4936").unwrap();
        let to_wgs84 = Transform::new("EPSG:4936", "EPSG:4978").unwrap();
        to_wgs84
            .convert_3d(to_etrs89.convert_3d(NYC_LON_LAT_H).unwrap())
            .unwrap()
    };
    assert!((from_compound.0 - via_etrs89_ecef.0).abs() < 1e-6);
    assert!((from_compound.1 - via_etrs89_ecef.1).abs() < 1e-6);
    assert!((from_compound.2 - via_etrs89_ecef.2).abs() < 1e-6);

    // True cross-datum: NAD27 geographic → ECEF matches the chained path.
    let nad27_direct = Transform::new("EPSG:4267", "EPSG:4978").unwrap();
    let nad27_chained = {
        let to_wgs84 = Transform::new("EPSG:4267", "EPSG:4326").unwrap();
        let to_ecef = Transform::new("EPSG:4326", "EPSG:4978").unwrap();
        let mid = to_wgs84.convert_3d((-90.0, 45.0, 250.0)).unwrap();
        to_ecef.convert_3d(mid).unwrap()
    };
    let nad27_ecef = nad27_direct.convert_3d((-90.0, 45.0, 250.0)).unwrap();
    assert!((nad27_ecef.0 - nad27_chained.0).abs() < 1e-6);
    assert!((nad27_ecef.1 - nad27_chained.1).abs() < 1e-6);
    assert!((nad27_ecef.2 - nad27_chained.2).abs() < 1e-6);
}

fn geographic_with_ellipsoidal_height(
    horizontal_epsg: u32,
    height_datum: Datum,
    unit: LinearUnit,
) -> CrsDef {
    CrsDef::Compound(Box::new(
        CompoundCrsDef::from_crs_def(
            0,
            proj_core::lookup_epsg(horizontal_epsg).unwrap(),
            VerticalCrsDef::ellipsoidal_height(0, height_datum, unit, ""),
            "",
        )
        .unwrap(),
    ))
}

#[test]
fn ellipsoidal_height_units_are_converted_at_ecef_boundary() {
    let source =
        geographic_with_ellipsoidal_height(4326, proj_core::datum::WGS84, LinearUnit::foot());
    let t = Transform::from_crs_defs(&source, &proj_core::lookup_epsg(4978).unwrap()).unwrap();

    let (x, _, _) = t.convert_3d((0.0, 0.0, 1.0)).unwrap();
    assert!((x - 6_378_137.304_8).abs() < 1e-6, "X = {x}");

    let (_, _, h) = t
        .inverse()
        .unwrap()
        .convert_3d((6_378_138.0, 0.0, 0.0))
        .unwrap();
    assert!((h - 1.0 / 0.3048).abs() < 1e-8, "h = {h}");
}

#[test]
fn mismatched_ellipsoidal_height_datum_to_ecef_is_rejected() {
    let source =
        geographic_with_ellipsoidal_height(4269, proj_core::datum::WGS84, LinearUnit::metre());
    let target = proj_core::lookup_epsg(4978).unwrap();
    let message = "ellipsoidal height datum must match the horizontal CRS datum";

    let forward = Transform::from_crs_defs(&source, &target).unwrap_err();
    assert!(forward.to_string().contains(message), "got {forward}");
    let inverse = Transform::from_crs_defs(&target, &source).unwrap_err();
    assert!(inverse.to_string().contains(message), "got {inverse}");
}

#[test]
fn gravity_related_height_to_ecef_is_rejected() {
    let err = Transform::new("EPSG:7415", "EPSG:4978").unwrap_err();
    assert!(
        err.to_string()
            .contains("explicit vertical CRS and a horizontal-only CRS"),
        "got {err}"
    );
}

#[test]
fn convert_2d_rejects_geocentric_endpoints() {
    let t = Transform::new("EPSG:4326", "EPSG:4978").unwrap();
    let err = t.convert((-74.006, 40.7128)).unwrap_err();
    assert!(err.to_string().contains("require convert_3d"), "got {err}");
    let err = t.convert_with_diagnostics((-74.006, 40.7128)).unwrap_err();
    assert!(err.to_string().contains("require convert_3d"), "got {err}");
}

#[test]
fn tuple3d_wgs84_to_web_mercator() {
    let t = Transform::new("EPSG:4326", "EPSG:3857").unwrap();
    let (x, y, z) = t.convert_3d((-74.0445, 40.6892, 15.5)).unwrap();

    assert!((x - (-8242596.0)).abs() < 1000.0);
    assert!((y - 4966606.0).abs() < 1000.0);
    assert!((z - 15.5).abs() < 1e-12);
}

#[test]
fn coord3d_roundtrip() {
    let fwd = Transform::new("EPSG:4326", "EPSG:3857").unwrap();
    let inv = Transform::new("EPSG:3857", "EPSG:4326").unwrap();

    let original = Coord3D::new(-74.0445, 40.6892, 25.0);
    let projected = fwd.convert_3d(original).unwrap();
    let roundtripped = inv.convert_3d(projected).unwrap();

    assert!((roundtripped.x - original.x).abs() < 1e-6);
    assert!((roundtripped.y - original.y).abs() < 1e-6);
    assert!((roundtripped.z - original.z).abs() < 1e-6);
}

#[test]
fn cross_datum_height_roundtrip() {
    let fwd = Transform::new("EPSG:4267", "EPSG:4326").unwrap();
    let inv = Transform::new("EPSG:4326", "EPSG:4267").unwrap();

    let original = (-90.0, 45.0, 250.0);
    let shifted = fwd.convert_3d(original).unwrap();
    let roundtripped = inv.convert_3d(shifted).unwrap();

    assert!((roundtripped.0 - original.0).abs() < 1e-6);
    assert!((roundtripped.1 - original.1).abs() < 1e-6);
    // Cartesian conversion around a Helmert step is accurate to well below a
    // micrometre, but does not promise picometre-exact geodetic heights.
    assert!((roundtripped.2 - original.2).abs() < 1e-8);
}

#[test]
fn helmert_backed_projected_transform_uses_source_height_for_xy() {
    let t = Transform::new("EPSG:4326", "EPSG:27700").unwrap();

    let ground = t.convert_3d((-0.1278, 51.5074, 0.0)).unwrap();
    let high = t.convert_3d((-0.1278, 51.5074, 10_000.0)).unwrap();
    let de = high.0 - ground.0;
    let dn = high.1 - ground.1;

    assert!(de.abs() > 0.1, "easting delta = {de}");
    assert!(dn.abs() > 0.05, "northing delta = {dn}");

    // OSGB36 Helmert is geographic-2D, so height is preserved.
    assert!(
        (ground.2 - 0.0).abs() < 1e-12,
        "ground height = {}",
        ground.2
    );
    assert!(
        (high.2 - 10_000.0).abs() < 1e-12,
        "high height = {}",
        high.2
    );
}

#[test]
fn wgs72_geocentric_to_wgs84_applies_full_3d_helmert() {
    // Both endpoints are geocentric, so Helmert runs in 3D.
    let ecef = Transform::new("EPSG:4322", "EPSG:4984")
        .unwrap()
        .convert_3d((0.0, 51.0, 100.0))
        .unwrap();
    let wgs84 = Transform::from_operation(CoordinateOperationId(1238), "EPSG:4984", "EPSG:4978")
        .unwrap()
        .convert_3d(ecef)
        .unwrap();
    let (_, _, h) = Transform::new("EPSG:4978", "EPSG:4326")
        .unwrap()
        .convert_3d(wgs84)
        .unwrap();
    assert!(
        (h - 100.0).abs() > 0.1,
        "3D Helmert should change height, got {h}"
    );
}

#[test]
fn rd_new_to_ecef_preserves_amersfoort_height() {
    // Amersfoort to WGS 84 (4) is geographic 2D; rotating height is ~43 m off libproj.
    let t = Transform::new("EPSG:28992", "EPSG:4978").unwrap();
    assert!(
        t.selected_operation().domain == OperationDomain::HorizontalOnly,
        "expected a horizontal-only operation, got {:?}",
        t.selected_operation().domain
    );

    let input = (155_000.0, 463_000.0, 43.0);
    let via_wgs84 = {
        let to_wgs84 = Transform::new("EPSG:28992", "EPSG:4326").unwrap();
        let to_ecef = Transform::new("EPSG:4326", "EPSG:4978").unwrap();
        let mid = to_wgs84.convert_3d(input).unwrap();
        assert!((mid.2 - 43.0).abs() < 1e-9, "geographic height = {}", mid.2);
        to_ecef.convert_3d(mid).unwrap()
    };
    let direct = t.convert_3d(input).unwrap();
    assert!((direct.0 - via_wgs84.0).abs() < 1e-6);
    assert!((direct.1 - via_wgs84.1).abs() < 1e-6);
    assert!((direct.2 - via_wgs84.2).abs() < 1e-6);
}

#[test]
fn batch_transform_3d_preserves_heights() {
    let t = Transform::new("EPSG:4326", "EPSG:3857").unwrap();
    let coords: Vec<(f64, f64, f64)> = (0..50)
        .map(|i| {
            (
                -74.0 + i as f64 * 0.01,
                40.0 + i as f64 * 0.01,
                i as f64 * 2.0,
            )
        })
        .collect();

    let results = t.convert_batch_3d(&coords).unwrap();
    assert_eq!(results.len(), 50);

    for (index, result) in results.iter().enumerate() {
        assert!(result.0 < 0.0);
        assert!((result.2 - index as f64 * 2.0).abs() < 1e-12);
    }
}

#[cfg(feature = "rayon")]
#[test]
fn parallel_batch_transform_3d_preserves_heights() {
    let t = Transform::new("EPSG:4326", "EPSG:3857").unwrap();
    let coords: Vec<(f64, f64, f64)> = (0..200)
        .map(|i| (-74.0 + i as f64 * 0.001, 40.0 + i as f64 * 0.001, i as f64))
        .collect();

    let results = t.convert_batch_parallel_3d(&coords).unwrap();
    assert_eq!(results.len(), 200);

    for (index, result) in results.iter().enumerate() {
        assert!((result.2 - index as f64).abs() < 1e-12);
    }
}
