#[path = "../src/viewer/alignment/dotplot/projection.rs"]
mod projection;
use projection::*;
fn gpu_axis(anchor: usize, camera: Camera, axis: usize) -> f64 {
    // Emulate the shader's f32 operations, then compare with independently
    // evaluated original-coordinate f64 base centers below.
    let local = (anchor as f32 - camera.integer[axis]) + (0.5f32 - camera.fraction[axis]);
    (local / camera.span[axis]) as f64
}
#[test]
fn one_base_zoom_preserves_high_coordinate_base_centers() {
    for fraction in [0.25, 0.75] {
        let origin = 9_999_998.0 + fraction;
        let c = camera(origin, origin, 1., 1.);
        for anchor in [9_999_998, 9_999_999] {
            let expected = anchor as f64 + 0.5 - origin;
            for axis in [0, 1] {
                assert!((gpu_axis(anchor, c, axis) - expected).abs() < 0.000001);
            }
        }
    }
    assert!(VERTEX_SHADER.contains("a_anchor.xy-u_origin_integer"));
}
#[test]
fn both_camera_axes_match_original_coordinates_across_zoom_levels() {
    for position in [0, 12345, 1_000_000, 9_999_998] {
        for fraction in [0., 0.125, 0.5, 0.875] {
            for span in [1., 4.5, 97., 1000.25] {
                let x = position as f64 + fraction;
                let y = (position / 2) as f64 + fraction;
                let c = camera(x, y, span, span * 2.);
                for delta in [0, 1] {
                    for (axis, anchor, origin, extent) in [
                        (0, position + delta, x, span),
                        (1, position / 2 + delta, y, span * 2.),
                    ] {
                        let expected = (anchor as f64 + 0.5 - origin) / extent;
                        assert!(
                            (gpu_axis(anchor, c, axis) - expected).abs() < 0.000001,
                            "{anchor} {origin} {extent}"
                        );
                    }
                }
            }
        }
    }
}
