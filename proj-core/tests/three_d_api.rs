use proj_core::{CompoundCrsDef, Coord3D, CrsDef, Datum, LinearUnit, Transform, VerticalCrsDef};

/// WGS 84 geodetic (lon/lat/h) → ECEF checkpoint for NYC from C PROJ /
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
fn wgs84_geographic_3d_to_ecef_roundtrip() {
    let t = Transform::new("EPSG:4979", "EPSG:4978").unwrap();
    let ecef = t.convert_3d(NYC_LON_LAT_H).unwrap();
    assert!((ecef.0 - NYC_ECEF_XYZ.0).abs() < 1e-4);
    assert!((ecef.1 - NYC_ECEF_XYZ.1).abs() < 1e-4);
    assert!((ecef.2 - NYC_ECEF_XYZ.2).abs() < 1e-4);

    let back = t.inverse().unwrap().convert_3d(ecef).unwrap();
    assert!((back.0 - NYC_LON_LAT_H.0).abs() < 1e-10);
    assert!((back.1 - NYC_LON_LAT_H.1).abs() < 1e-10);
    assert!((back.2 - NYC_LON_LAT_H.2).abs() < 1e-6);
}

#[test]
fn cross_datum_ellipsoidal_3d_to_ecef() {
    // ETRS89 3D (ellipsoidal height) → WGS 84 ECEF must be allowed: height is
    // consumed by cart framing after the selected horizontal datum operation.
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

    // Same horizontal datum without a compound vertical CRS should match:
    // convert_3d treats z as ellipsoidal height into ECEF.
    let geographic = Transform::new("EPSG:4258", "EPSG:4978").unwrap();
    let from_compound = compound.convert_3d(NYC_LON_LAT_H).unwrap();
    let from_geographic = geographic.convert_3d(NYC_LON_LAT_H).unwrap();
    assert!((from_compound.0 - from_geographic.0).abs() < 1e-9);
    assert!((from_compound.1 - from_geographic.1).abs() < 1e-9);
    assert!((from_compound.2 - from_geographic.2).abs() < 1e-9);

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

#[test]
fn projected_to_ecef_matches_geographic_path() {
    let to_utm = Transform::new("EPSG:4326", "EPSG:32618").unwrap();
    let utm = to_utm.convert_3d(NYC_LON_LAT_H).unwrap();

    let utm_to_ecef = Transform::new("EPSG:32618", "EPSG:4978").unwrap();
    let ecef = utm_to_ecef.convert_3d(utm).unwrap();
    assert!((ecef.0 - NYC_ECEF_XYZ.0).abs() < 1e-3, "x = {}", ecef.0);
    assert!((ecef.1 - NYC_ECEF_XYZ.1).abs() < 1e-3, "y = {}", ecef.1);
    assert!((ecef.2 - NYC_ECEF_XYZ.2).abs() < 1e-3, "z = {}", ecef.2);

    let back = utm_to_ecef.inverse().unwrap().convert_3d(ecef).unwrap();
    assert!((back.0 - utm.0).abs() < 1e-6);
    assert!((back.1 - utm.1).abs() < 1e-6);
    assert!((back.2 - utm.2).abs() < 1e-6);
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
    let to_ecef = Transform::new("EPSG:4326", "EPSG:4978").unwrap();
    let err = to_ecef.convert((-74.006, 40.7128)).unwrap_err();
    assert!(err.to_string().contains("require convert_3d"), "got {err}");

    let from_ecef = Transform::new("EPSG:4978", "EPSG:4326").unwrap();
    let err = from_ecef
        .convert((NYC_ECEF_XYZ.0, NYC_ECEF_XYZ.1))
        .unwrap_err();
    assert!(err.to_string().contains("require convert_3d"), "got {err}");

    let err = to_ecef
        .convert_with_diagnostics((-74.006, 40.7128))
        .unwrap_err();
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

    // Ellipsoidal height is rebased onto the target datum's ellipsoid: the
    // WGS84→OSGB36(Airy) separation near London is about -46 m, and height
    // differences are preserved to within the datum shift's height gradient.
    assert!(
        (-60.0..=-30.0).contains(&ground.2),
        "ground height = {}",
        ground.2
    );
    assert!(
        (high.2 - ground.2 - 10_000.0).abs() < 20.0,
        "height delta = {}",
        high.2 - ground.2
    );
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
