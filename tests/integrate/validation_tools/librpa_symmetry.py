"""Small matrix checks for LibRPA's row-vector symmetry metadata."""


def reciprocal_rotation(direct_rotation):
    """Return the inverse transpose for a direct-space row-vector operation."""
    if len(direct_rotation) != 9:
        raise ValueError("Symmetry rotation must have nine elements")
    a, b, c, d, e, f, g, h, i = map(int, direct_rotation)
    determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    if abs(determinant) != 1:
        raise ValueError("Symmetry rotation must be unimodular")
    return [
        determinant * (e * i - f * h),
        determinant * (f * g - d * i),
        determinant * (d * h - e * g),
        determinant * (c * h - b * i),
        determinant * (a * i - c * g),
        determinant * (b * g - a * h),
        determinant * (b * f - c * e),
        determinant * (c * d - a * f),
        determinant * (a * e - b * d),
    ]


def preserves_bvk_mesh(direct_rotation, mesh, offset=(0.0, 0.0, 0.0), antiunitary=False):
    """Check that a spatial or antiunitary operation maps an MP mesh to itself."""
    reciprocal = reciprocal_rotation(direct_rotation)
    sign = -1 if antiunitary else 1
    for source in range(3):
        for target in range(3):
            numerator = sign * reciprocal[3 * source + target] * mesh[target]
            if numerator % mesh[source] != 0:
                return False

    origin = [offset[index] / mesh[index] for index in range(3)]
    for target in range(3):
        mapped = sign * sum(origin[source] * reciprocal[3 * source + target] for source in range(3))
        shift = (mapped - origin[target]) * mesh[target]
        if abs(shift - round(shift)) > 1e-12:
            return False
    return True
