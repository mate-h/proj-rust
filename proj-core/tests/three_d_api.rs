use proj_core::{Coord3D, Transform};

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
