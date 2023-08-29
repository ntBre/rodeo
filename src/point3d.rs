use std::ops::Sub;

#[derive(Clone, Debug, Default)]
#[allow(unused)]
pub(crate) struct Point3D {
    pub(crate) x: f64,
    pub(crate) y: f64,
    pub(crate) z: f64,
}

impl Point3D {
    pub(crate) fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    #[inline]
    pub(crate) const fn zero() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    #[inline]
    pub(crate) fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub(crate) fn cross(&self, other: &Self) -> Self {
        let Self {
            x: a1,
            y: a2,
            z: a3,
        } = self;
        let Self {
            x: b1,
            y: b2,
            z: b3,
        } = other;
        Self {
            x: a2 * b3 - a3 * b2,
            y: a3 * b1 - a1 * b3,
            z: a1 * b2 - a2 * b1,
        }
    }
}

impl Sub<&Point3D> for Point3D {
    type Output = Self;

    fn sub(mut self, rhs: &Point3D) -> Self::Output {
        self.x -= rhs.x;
        self.y -= rhs.y;
        self.z -= rhs.z;
        self
    }
}
