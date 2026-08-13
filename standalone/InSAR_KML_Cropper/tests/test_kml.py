from insar_kml_cropper.kml import expand_bounds


def test_expand_bounds_twenty_percent_per_side():
    assert expand_bounds((10.0, 20.0, 30.0, 40.0), 0.2) == (8.0, 22.0, 28.0, 42.0)
