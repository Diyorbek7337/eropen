
class Sphere {
	constructor(e = new Vector3, t = -1) {
		this.center = e, this.radius = t
	}
	set(e, t) {
		return this.center.copy(e), this.radius = t, this
	}
	setFromPoints(e, t) {
		const i = this.center;
		void 0 !== t ? i.copy(t) : _box$1.setFromPoints(e).getCenter(i);
		let n = 0;
		for (let t = 0, r = e.length; t < r; t++) n = Math.max(n, i.distanceToSquared(e[t]));
		return this.radius = Math.sqrt(n), this
	}
	copy(e) {
		return this.center.copy(e.center), this.radius = e.radius, this
	}
	isEmpty() {
		return this.radius < 0
	}
	makeEmpty() {
		return this.center.set(0, 0, 0), this.radius = -1, this
	}
	containsPoint(e) {
		return e.distanceToSquared(this.center) <= this.radius * this.radius
	}
	distanceToPoint(e) {
		return e.distanceTo(this.center) - this.radius
	}
	intersectsSphere(e) {
		const t = this.radius + e.radius;
		return e.center.distanceToSquared(this.center) <= t * t
	}
	intersectsBox(e) {
		return e.intersectsSphere(this)
	}
	intersectsPlane(e) {
		return Math.abs(e.distanceToPoint(this.center)) <= this.radius
	}
	clampPoint(e, t) {
		const i = this.center.distanceToSquared(e);
		return void 0 === t && (console.warn("THREE.Sphere: .clampPoint() target is now required"), t = new Vector3), t.copy(e), i > this.radius * this.radius && (t.sub(this.center).normalize(), t.multiplyScalar(this.radius).add(this.center)), t
	}
	getBoundingBox(e) {
		return void 0 === e && (console.warn("THREE.Sphere: .getBoundingBox() target is now required"), e = new Box3), this.isEmpty() ? (e.makeEmpty(), e) : (e.set(this.center, this.center), e.expandByScalar(this.radius), e)
	}
	applyMatrix4(e) {
		return this.center.applyMatrix4(e), this.radius = this.radius * e.getMaxScaleOnAxis(), this
	}
	translate(e) {
		return this.center.add(e), this
	}
	equals(e) {
		return e.center.equals(this.center) && e.radius === this.radius
	}
	clone() {
		return (new this.constructor).copy(this)
	}
}
const _vector$2 = new Vector3,
	_segCenter = new Vector3,
	_segDir = new Vector3,
	_diff = new Vector3,
	_edge1 = new Vector3,
	_edge2 = new Vector3,
	_normal = new Vector3;
class Ray {
	constructor(e = new Vector3, t = new Vector3(0, 0, -1)) {
		this.origin = e, this.direction = t
	}
	set(e, t) {
		return this.origin.copy(e), this.direction.copy(t), this
	}
	copy(e) {
		return this.origin.copy(e.origin), this.direction.copy(e.direction), this
	}
	at(e, t) {
		return void 0 === t && (console.warn("THREE.Ray: .at() target is now required"), t = new Vector3), t.copy(this.direction).multiplyScalar(e).add(this.origin)
	}
	lookAt(e) {
		return this.direction.copy(e).sub(this.origin).normalize(), this
	}
	recast(e) {
		return this.origin.copy(this.at(e, _vector$2)), this
	}
	closestPointToPoint(e, t) {
		void 0 === t && (console.warn("THREE.Ray: .closestPointToPoint() target is now required"), t = new Vector3), t.subVectors(e, this.origin);
		const i = t.dot(this.direction);
		return i < 0 ? t.copy(this.origin) : t.copy(this.direction).multiplyScalar(i).add(this.origin)
	}
	distanceToPoint(e) {
		return Math.sqrt(this.distanceSqToPoint(e))
	}
	distanceSqToPoint(e) {
		const t = _vector$2.subVectors(e, this.origin).dot(this.direction);
		return t < 0 ? this.origin.distanceToSquared(e) : (_vector$2.copy(this.direction).multiplyScalar(t).add(this.origin), _vector$2.distanceToSquared(e))
	}
	distanceSqToSegment(e, t, i, n) {
		_segCenter.copy(e).add(t).multiplyScalar(.5), _segDir.copy(t).sub(e).normalize(), _diff.copy(this.origin).sub(_segCenter);
		const r = .5 * e.distanceTo(t),
			a = -this.direction.dot(_segDir),
			s = _diff.dot(this.direction),
			o = -_diff.dot(_segDir),
			l = _diff.lengthSq(),
			c = Math.abs(1 - a * a);
		let h, u, d, p;
		if (c > 0)
			if (h = a * o - s, u = a * s - o, p = r * c, h >= 0)
				if (u >= -p)
					if (u <= p) {
						const e = 1 / c;
						h *= e, u *= e, d = h * (h + a * u + 2 * s) + u * (a * h + u + 2 * o) + l
					} else u = r, h = Math.max(0, -(a * u + s)), d = -h * h + u * (u + 2 * o) + l;
		else u = -r, h = Math.max(0, -(a * u + s)), d = -h * h + u * (u + 2 * o) + l;
		else u <= -p ? (h = Math.max(0, -(-a * r + s)), u = h > 0 ? -r : Math.min(Math.max(-r, -o), r), d = -h * h + u * (u + 2 * o) + l) : u <= p ? (h = 0, u = Math.min(Math.max(-r, -o), r), d = u * (u + 2 * o) + l) : (h = Math.max(0, -(a * r + s)), u = h > 0 ? r : Math.min(Math.max(-r, -o), r), d = -h * h + u * (u + 2 * o) + l);
		else u = a > 0 ? -r : r, h = Math.max(0, -(a * u + s)), d = -h * h + u * (u + 2 * o) + l;
		return i && i.copy(this.direction).multiplyScalar(h).add(this.origin), n && n.copy(_segDir).multiplyScalar(u).add(_segCenter), d
	}
	intersectSphere(e, t) {
		_vector$2.subVectors(e.center, this.origin);
		const i = _vector$2.dot(this.direction),
			n = _vector$2.dot(_vector$2) - i * i,
			r = e.radius * e.radius;
		if (n > r) return null;
		const a = Math.sqrt(r - n),
			s = i - a,
			o = i + a;
		return s < 0 && o < 0 ? null : s < 0 ? this.at(o, t) : this.at(s, t)
	}
	intersectsSphere(e) {
		return this.distanceSqToPoint(e.center) <= e.radius * e.radius
	}
	distanceToPlane(e) {
		const t = e.normal.dot(this.direction);
		if (0 === t) return 0 === e.distanceToPoint(this.origin) ? 0 : null;
		const i = -(this.origin.dot(e.normal) + e.constant) / t;
		return i >= 0 ? i : null
	}
	intersectPlane(e, t) {
		const i = this.distanceToPlane(e);
		return null === i ? null : this.at(i, t)
	}
	intersectsPlane(e) {
		const t = e.distanceToPoint(this.origin);
		if (0 === t) return !0;
		return e.normal.dot(this.direction) * t < 0
	}
	intersectBox(e, t) {
		let i, n, r, a, s, o;
		const l = 1 / this.direction.x,
			c = 1 / this.direction.y,
			h = 1 / this.direction.z,
			u = this.origin;
		return l >= 0 ? (i = (e.min.x - u.x) * l, n = (e.max.x - u.x) * l) : (i = (e.max.x - u.x) * l, n = (e.min.x - u.x) * l), c >= 0 ? (r = (e.min.y - u.y) * c, a = (e.max.y - u.y) * c) : (r = (e.max.y - u.y) * c, a = (e.min.y - u.y) * c), i > a || r > n ? null : ((r > i || i != i) && (i = r), (a < n || n != n) && (n = a), h >= 0 ? (s = (e.min.z - u.z) * h, o = (e.max.z - u.z) * h) : (s = (e.max.z - u.z) * h, o = (e.min.z - u.z) * h), i > o || s > n ? null : ((s > i || i != i) && (i = s), (o < n || n != n) && (n = o), n < 0 ? null : this.at(i >= 0 ? i : n, t)))
	}
	intersectsBox(e) {
		return null !== this.intersectBox(e, _vector$2)
	}
	intersectTriangle(e, t, i, n, r) {
		_edge1.subVectors(t, e), _edge2.subVectors(i, e), _normal.crossVectors(_edge1, _edge2);
		let a, s = this.direction.dot(_normal);
		if (s > 0) {
			if (n) return null;
			a = 1
		} else {
			if (!(s < 0)) return null;
			a = -1, s = -s
		}
		_diff.subVectors(this.origin, e);
		const o = a * this.direction.dot(_edge2.crossVectors(_diff, _edge2));
		if (o < 0) return null;
		const l = a * this.direction.dot(_edge1.cross(_diff));
		if (l < 0) return null;
		if (o + l > s) return null;
		const c = -a * _diff.dot(_normal);
		return c < 0 ? null : this.at(c / s, r)
	}
	applyMatrix4(e) {
		return this.origin.applyMatrix4(e), this.direction.transformDirection(e), this
	}
	equals(e) {
		return e.origin.equals(this.origin) && e.direction.equals(this.direction)
	}
	clone() {
		return (new this.constructor).copy(this)
	}
}
class Matrix4 {
	constructor() {
		this.elements = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1], arguments.length > 0 && console.error("THREE.Matrix4: the constructor no longer reads arguments. use .set() instead.")
	}
	set(e, t, i, n, r, a, s, o, l, c, h, u, d, p, m, A) {
		const g = this.elements;
		return g[0] = e, g[4] = t, g[8] = i, g[12] = n, g[1] = r, g[5] = a, g[9] = s, g[13] = o, g[2] = l, g[6] = c, g[10] = h, g[14] = u, g[3] = d, g[7] = p, g[11] = m, g[15] = A, this
	}
	identity() {
		return this.set(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), this
	}
	clone() {
		return (new Matrix4).fromArray(this.elements)
	}
	copy(e) {
		const t = this.elements,
			i = e.elements;
		return t[0] = i[0], t[1] = i[1], t[2] = i[2], t[3] = i[3], t[4] = i[4], t[5] = i[5], t[6] = i[6], t[7] = i[7], t[8] = i[8], t[9] = i[9], t[10] = i[10], t[11] = i[11], t[12] = i[12], t[13] = i[13], t[14] = i[14], t[15] = i[15], this
	}
	copyPosition(e) {
		const t = this.elements,
			i = e.elements;
		return t[12] = i[12], t[13] = i[13], t[14] = i[14], this
	}
	setFromMatrix3(e) {
		const t = e.elements;
		return this.set(t[0], t[3], t[6], 0, t[1], t[4], t[7], 0, t[2], t[5], t[8], 0, 0, 0, 0, 1), this
	}
	extractBasis(e, t, i) {
		return e.setFromMatrixColumn(this, 0), t.setFromMatrixColumn(this, 1), i.setFromMatrixColumn(this, 2), this
	}
	makeBasis(e, t, i) {
		return this.set(e.x, t.x, i.x, 0, e.y, t.y, i.y, 0, e.z, t.z, i.z, 0, 0, 0, 0, 1), this
	}
	extractRotation(e) {
		const t = this.elements,
			i = e.elements,
			n = 1 / _v1$1.setFromMatrixColumn(e, 0).length(),
			r = 1 / _v1$1.setFromMatrixColumn(e, 1).length(),
			a = 1 / _v1$1.setFromMatrixColumn(e, 2).length();
		return t[0] = i[0] * n, t[1] = i[1] * n, t[2] = i[2] * n, t[3] = 0, t[4] = i[4] * r, t[5] = i[5] * r, t[6] = i[6] * r, t[7] = 0, t[8] = i[8] * a, t[9] = i[9] * a, t[10] = i[10] * a, t[11] = 0, t[12] = 0, t[13] = 0, t[14] = 0, t[15] = 1, this
	}
	makeRotationFromEuler(e) {
		e && e.isEuler || console.error("THREE.Matrix4: .makeRotationFromEuler() now expects a Euler rotation rather than a Vector3 and order.");
		const t = this.elements,
			i = e.x,
			n = e.y,
			r = e.z,
			a = Math.cos(i),
			s = Math.sin(i),
			o = Math.cos(n),
			l = Math.sin(n),
			c = Math.cos(r),
			h = Math.sin(r);
		if ("XYZ" === e.order) {
			const e = a * c,
				i = a * h,
				n = s * c,
				r = s * h;
			t[0] = o * c, t[4] = -o * h, t[8] = l, t[1] = i + n * l, t[5] = e - r * l, t[9] = -s * o, t[2] = r - e * l, t[6] = n + i * l, t[10] = a * o
		} else if ("YXZ" === e.order) {
			const e = o * c,
				i = o * h,
				n = l * c,
				r = l * h;
			t[0] = e + r * s, t[4] = n * s - i, t[8] = a * l, t[1] = a * h, t[5] = a * c, t[9] = -s, t[2] = i * s - n, t[6] = r + e * s, t[10] = a * o
		} else if ("ZXY" === e.order) {
			const e = o * c,
				i = o * h,
				n = l * c,
				r = l * h;
			t[0] = e - r * s, t[4] = -a * h, t[8] = n + i * s, t[1] = i + n * s, t[5] = a * c, t[9] = r - e * s, t[2] = -a * l, t[6] = s, t[10] = a * o
		} else if ("ZYX" === e.order) {
			const e = a * c,
				i = a * h,
				n = s * c,
				r = s * h;
			t[0] = o * c, t[4] = n * l - i, t[8] = e * l + r, t[1] = o * h, t[5] = r * l + e, t[9] = i * l - n, t[2] = -l, t[6] = s * o, t[10] = a * o
		} else if ("YZX" === e.order) {
			const e = a * o,
				i = a * l,
				n = s * o,
				r = s * l;
			t[0] = o * c, t[4] = r - e * h, t[8] = n * h + i, t[1] = h, t[5] = a * c, t[9] = -s * c, t[2] = -l * c, t[6] = i * h + n, t[10] = e - r * h
		} else if ("XZY" === e.order) {
			const e = a * o,
				i = a * l,
				n = s * o,
				r = s * l;
			t[0] = o * c, t[4] = -h, t[8] = l * c, t[1] = e * h + r, t[5] = a * c, t[9] = i * h - n, t[2] = n * h - i, t[6] = s * c, t[10] = r * h + e
		}
		return t[3] = 0, t[7] = 0, t[11] = 0, t[12] = 0, t[13] = 0, t[14] = 0, t[15] = 1, this
	}
	makeRotationFromQuaternion(e) {
		return this.compose(_zero, e, _one)
	}
	lookAt(e, t, i) {
		const n = this.elements;
		return _z.subVectors(e, t), 0 === _z.lengthSq() && (_z.z = 1), _z.normalize(), _x.crossVectors(i, _z), 0 === _x.lengthSq() && (1 === Math.abs(i.z) ? _z.x += 1e-4 : _z.z += 1e-4, _z.normalize(), _x.crossVectors(i, _z)), _x.normalize(), _y.crossVectors(_z, _x), n[0] = _x.x, n[4] = _y.x, n[8] = _z.x, n[1] = _x.y, n[5] = _y.y, n[9] = _z.y, n[2] = _x.z, n[6] = _y.z, n[10] = _z.z, this
	}
	multiply(e, t) {
		return void 0 !== t ? (console.warn("THREE.Matrix4: .multiply() now only accepts one argument. Use .multiplyMatrices( a, b ) instead."), this.multiplyMatrices(e, t)) : this.multiplyMatrices(this, e)
	}
	premultiply(e) {
		return this.multiplyMatrices(e, this)
	}
	multiplyMatrices(e, t) {
		const i = e.elements,
			n = t.elements,
			r = this.elements,
			a = i[0],
			s = i[4],
			o = i[8],
			l = i[12],
			c = i[1],
			h = i[5],
			u = i[9],
			d = i[13],
			p = i[2],
			m = i[6],
			A = i[10],
			g = i[14],
			f = i[3],
			v = i[7],
			y = i[11],
			E = i[15],
			_ = n[0],
			b = n[4],
			x = n[8],
			w = n[12],
			C = n[1],
			S = n[5],
			I = n[9],
			M = n[13],
			T = n[2],
			B = n[6],
			L = n[10],
			R = n[14],
			D = n[3],
			P = n[7],
			Q = n[11],
			F = n[15];
		return r[0] = a * _ + s * C + o * T + l * D, r[4] = a * b + s * S + o * B + l * P, r[8] = a * x + s * I + o * L + l * Q, r[12] = a * w + s * M + o * R + l * F, r[1] = c * _ + h * C + u * T + d * D, r[5] = c * b + h * S + u * B + d * P, r[9] = c * x + h * I + u * L + d * Q, r[13] = c * w + h * M + u * R + d * F, r[2] = p * _ + m * C + A * T + g * D, r[6] = p * b + m * S + A * B + g * P, r[10] = p * x + m * I + A * L + g * Q, r[14] = p * w + m * M + A * R + g * F, r[3] = f * _ + v * C + y * T + E * D, r[7] = f * b + v * S + y * B + E * P, r[11] = f * x + v * I + y * L + E * Q, r[15] = f * w + v * M + y * R + E * F, this
	}
	multiplyScalar(e) {
		const t = this.elements;
		return t[0] *= e, t[4] *= e, t[8] *= e, t[12] *= e, t[1] *= e, t[5] *= e, t[9] *= e, t[13] *= e, t[2] *= e, t[6] *= e, t[10] *= e, t[14] *= e, t[3] *= e, t[7] *= e, t[11] *= e, t[15] *= e, this
	}
	determinant() {
		const e = this.elements,
			t = e[0],
			i = e[4],
			n = e[8],
			r = e[12],
			a = e[1],
			s = e[5],
			o = e[9],
			l = e[13],
			c = e[2],
			h = e[6],
			u = e[10],
			d = e[14];
		return e[3] * (+r * o * h - n * l * h - r * s * u + i * l * u + n * s * d - i * o * d) + e[7] * (+t * o * d - t * l * u + r * a * u - n * a * d + n * l * c - r * o * c) + e[11] * (+t * l * h - t * s * d - r * a * h + i * a * d + r * s * c - i * l * c) + e[15] * (-n * s * c - t * o * h + t * s * u + n * a * h - i * a * u + i * o * c)
	}
	transpose() {
		const e = this.elements;
		let t;
		return t = e[1], e[1] = e[4], e[4] = t, t = e[2], e[2] = e[8], e[8] = t, t = e[6], e[6] = e[9], e[9] = t, t = e[3], e[3] = e[12], e[12] = t, t = e[7], e[7] = e[13], e[13] = t, t = e[11], e[11] = e[14], e[14] = t, this
	}
	setPosition(e, t, i) {
		const n = this.elements;
		return e.isVector3 ? (n[12] = e.x, n[13] = e.y, n[14] = e.z) : (n[12] = e, n[13] = t, n[14] = i), this
	}
	invert() {
		const e = this.elements,
			t = e[0],
			i = e[1],
			n = e[2],
			r = e[3],
			a = e[4],
			s = e[5],
			o = e[6],
			l = e[7],
			c = e[8],
			h = e[9],
			u = e[10],
			d = e[11],
			p = e[12],
			m = e[13],
			A = e[14],
			g = e[15],
			f = h * A * l - m * u * l + m * o * d - s * A * d - h * o * g + s * u * g,
			v = p * u * l - c * A * l - p * o * d + a * A * d + c * o * g - a * u * g,
			y = c * m * l - p * h * l + p * s * d - a * m * d - c * s * g + a * h * g,
			E = p * h * o - c * m * o - p * s * u + a * m * u + c * s * A - a * h * A,
			_ = t * f + i * v + n * y + r * E;
		if (0 === _) return this.set(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
		const b = 1 / _;
		return e[0] = f * b, e[1] = (m * u * r - h * A * r - m * n * d + i * A * d + h * n * g - i * u * g) * b, e[2] = (s * A * r - m * o * r + m * n * l - i * A * l - s * n * g + i * o * g) * b, e[3] = (h * o * r - s * u * r - h * n * l + i * u * l + s * n * d - i * o * d) * b, e[4] = v * b, e[5] = (c * A * r - p * u * r + p * n * d - t * A * d - c * n * g + t * u * g) * b, e[6] = (p * o * r - a * A * r - p * n * l + t * A * l + a * n * g - t * o * g) * b, e[7] = (a * u * r - c * o * r + c * n * l - t * u * l - a * n * d + t * o * d) * b, e[8] = y * b, e[9] = (p * h * r - c * m * r - p * i * d + t * m * d + c * i * g - t * h * g) * b, e[10] = (a * m * r - p * s * r + p * i * l - t * m * l - a * i * g + t * s * g) * b, e[11] = (c * s * r - a * h * r - c * i * l + t * h * l + a * i * d - t * s * d) * b, e[12] = E * b, e[13] = (c * m * n - p * h * n + p * i * u - t * m * u - c * i * A + t * h * A) * b, e[14] = (p * s * n - a * m * n - p * i * o + t * m * o + a * i * A - t * s * A) * b, e[15] = (a * h * n - c * s * n + c * i * o - t * h * o - a * i * u + t * s * u) * b, this
	}
	scale(e) {
		const t = this.elements,
			i = e.x,
			n = e.y,
			r = e.z;
		return t[0] *= i, t[4] *= n, t[8] *= r, t[1] *= i, t[5] *= n, t[9] *= r, t[2] *= i, t[6] *= n, t[10] *= r, t[3] *= i, t[7] *= n, t[11] *= r, this
	}
	getMaxScaleOnAxis() {
		const e = this.elements,
			t = e[0] * e[0] + e[1] * e[1] + e[2] * e[2],
			i = e[4] * e[4] + e[5] * e[5] + e[6] * e[6],
			n = e[8] * e[8] + e[9] * e[9] + e[10] * e[10];
		return Math.sqrt(Math.max(t, i, n))
	}
	makeTranslation(e, t, i) {
		return this.set(1, 0, 0, e, 0, 1, 0, t, 0, 0, 1, i, 0, 0, 0, 1), this
	}
	makeRotationX(e) {
		const t = Math.cos(e),
			i = Math.sin(e);
		return this.set(1, 0, 0, 0, 0, t, -i, 0, 0, i, t, 0, 0, 0, 0, 1), this
	}
	makeRotationY(e) {
		const t = Math.cos(e),
			i = Math.sin(e);
		return this.set(t, 0, i, 0, 0, 1, 0, 0, -i, 0, t, 0, 0, 0, 0, 1), this
	}
	makeRotationZ(e) {
		const t = Math.cos(e),
			i = Math.sin(e);
		return this.set(t, -i, 0, 0, i, t, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), this
	}
	makeRotationAxis(e, t) {
		const i = Math.cos(t),
			n = Math.sin(t),
			r = 1 - i,
			a = e.x,
			s = e.y,
			o = e.z,
			l = r * a,
			c = r * s;
		return this.set(l * a + i, l * s - n * o, l * o + n * s, 0, l * s + n * o, c * s + i, c * o - n * a, 0, l * o - n * s, c * o + n * a, r * o * o + i, 0, 0, 0, 0, 1), this
	}
	makeScale(e, t, i) {
		return this.set(e, 0, 0, 0, 0, t, 0, 0, 0, 0, i, 0, 0, 0, 0, 1), this
	}
	makeShear(e, t, i) {
		return this.set(1, t, i, 0, e, 1, i, 0, e, t, 1, 0, 0, 0, 0, 1), this
	}
	compose(e, t, i) {
		const n = this.elements,
			r = t._x,
			a = t._y,
			s = t._z,
			o = t._w,
			l = r + r,
			c = a + a,
			h = s + s,
			u = r * l,
			d = r * c,
			p = r * h,
			m = a * c,
			A = a * h,
			g = s * h,
			f = o * l,
			v = o * c,
			y = o * h,
			E = i.x,
			_ = i.y,
			b = i.z;
		return n[0] = (1 - (m + g)) * E, n[1] = (d + y) * E, n[2] = (p - v) * E, n[3] = 0, n[4] = (d - y) * _, n[5] = (1 - (u + g)) * _, n[6] = (A + f) * _, n[7] = 0, n[8] = (p + v) * b, n[9] = (A - f) * b, n[10] = (1 - (u + m)) * b, n[11] = 0, n[12] = e.x, n[13] = e.y, n[14] = e.z, n[15] = 1, this
	}
	decompose(e, t, i) {
		const n = this.elements;
		let r = _v1$1.set(n[0], n[1], n[2]).length();
		const a = _v1$1.set(n[4], n[5], n[6]).length(),
			s = _v1$1.set(n[8], n[9], n[10]).length();
		this.determinant() < 0 && (r = -r), e.x = n[12], e.y = n[13], e.z = n[14], _m1.copy(this);
		const o = 1 / r,
			l = 1 / a,
			c = 1 / s;
		return _m1.elements[0] *= o, _m1.elements[1] *= o, _m1.elements[2] *= o, _m1.elements[4] *= l, _m1.elements[5] *= l, _m1.elements[6] *= l, _m1.elements[8] *= c, _m1.elements[9] *= c, _m1.elements[10] *= c, t.setFromRotationMatrix(_m1), i.x = r, i.y = a, i.z = s, this
	}
	makePerspective(e, t, i, n, r, a) {
		void 0 === a && console.warn("THREE.Matrix4: .makePerspective() has been redefined and has a new signature. Please check the docs.");
		const s = this.elements,
			o = 2 * r / (t - e),
			l = 2 * r / (i - n),
			c = (t + e) / (t - e),
			h = (i + n) / (i - n),
			u = -(a + r) / (a - r),
			d = -2 * a * r / (a - r);
		return s[0] = o, s[4] = 0, s[8] = c, s[12] = 0, s[1] = 0, s[5] = l, s[9] = h, s[13] = 0, s[2] = 0, s[6] = 0, s[10] = u, s[14] = d, s[3] = 0, s[7] = 0, s[11] = -1, s[15] = 0, this
	}
	makeOrthographic(e, t, i, n, r, a) {
		const s = this.elements,
			o = 1 / (t - e),
			l = 1 / (i - n),
			c = 1 / (a - r),
			h = (t + e) * o,
			u = (i + n) * l,
			d = (a + r) * c;
		return s[0] = 2 * o, s[4] = 0, s[8] = 0, s[12] = -h, s[1] = 0, s[5] = 2 * l, s[9] = 0, s[13] = -u, s[2] = 0, s[6] = 0, s[10] = -2 * c, s[14] = -d, s[3] = 0, s[7] = 0, s[11] = 0, s[15] = 1, this
	}
	equals(e) {
		const t = this.elements,
			i = e.elements;
		for (let e = 0; e < 16; e++)
			if (t[e] !== i[e]) return !1;
		return !0
	}
	fromArray(e, t = 0) {
		for (let i = 0; i < 16; i++) this.elements[i] = e[i + t];
		return this
	}
	toArray(e = [], t = 0) {
		const i = this.elements;
		return e[t] = i[0], e[t + 1] = i[1], e[t + 2] = i[2], e[t + 3] = i[3], e[t + 4] = i[4], e[t + 5] = i[5], e[t + 6] = i[6], e[t + 7] = i[7], e[t + 8] = i[8], e[t + 9] = i[9], e[t + 10] = i[10], e[t + 11] = i[11], e[t + 12] = i[12], e[t + 13] = i[13], e[t + 14] = i[14], e[t + 15] = i[15], e
	}
}
Matrix4.prototype.isMatrix4 = !0;
const _v1$1 = new Vector3,
	_m1 = new Matrix4,
	_zero = new Vector3(0, 0, 0),
	_one = new Vector3(1, 1, 1),
	_x = new Vector3,
	_y = new Vector3,
	_z = new Vector3,
	_matrix = new Matrix4,
	_quaternion$1 = new Quaternion;
class Euler {
	constructor(e = 0, t = 0, i = 0, n = Euler.DefaultOrder) {
		this._x = e, this._y = t, this._z = i, this._order = n
	}
	get x() {
		return this._x
	}
	set x(e) {
		this._x = e, this._onChangeCallback()
	}
	get y() {
		return this._y
	}
	set y(e) {
		this._y = e, this._onChangeCallback()
	}
	get z() {
		return this._z
	}
	set z(e) {
		this._z = e, this._onChangeCallback()
	}
	get order() {
		return this._order
	}
	set order(e) {
		this._order = e, this._onChangeCallback()
	}
	set(e, t, i, n) {
		return this._x = e, this._y = t, this._z = i, this._order = n || this._order, this._onChangeCallback(), this
	}
	clone() {
		return new this.constructor(this._x, this._y, this._z, this._order)
	}
	copy(e) {
		return this._x = e._x, this._y = e._y, this._z = e._z, this._order = e._order, this._onChangeCallback(), this
	}
	setFromRotationMatrix(e, t, i) {
		const n = MathUtils.clamp,
			r = e.elements,
			a = r[0],
			s = r[4],
			o = r[8],
			l = r[1],
			c = r[5],
			h = r[9],
			u = r[2],
			d = r[6],
			p = r[10];
		switch (t = t || this._order) {
			case "XYZ":
				this._y = Math.asin(n(o, -1, 1)), Math.abs(o) < .9999999 ? (this._x = Math.atan2(-h, p), this._z = Math.atan2(-s, a)) : (this._x = Math.atan2(d, c), this._z = 0);
				break;
			case "YXZ":
				this._x = Math.asin(-n(h, -1, 1)), Math.abs(h) < .9999999 ? (this._y = Math.atan2(o, p), this._z = Math.atan2(l, c)) : (this._y = Math.atan2(-u, a), this._z = 0);
				break;
			case "ZXY":
				this._x = Math.asin(n(d, -1, 1)), Math.abs(d) < .9999999 ? (this._y = Math.atan2(-u, p), this._z = Math.atan2(-s, c)) : (this._y = 0, this._z = Math.atan2(l, a));
				break;
			case "ZYX":
				this._y = Math.asin(-n(u, -1, 1)), Math.abs(u) < .9999999 ? (this._x = Math.atan2(d, p), this._z = Math.atan2(l, a)) : (this._x = 0, this._z = Math.atan2(-s, c));
				break;
			case "YZX":
				this._z = Math.asin(n(l, -1, 1)), Math.abs(l) < .9999999 ? (this._x = Math.atan2(-h, c), this._y = Math.atan2(-u, a)) : (this._x = 0, this._y = Math.atan2(o, p));
				break;
			case "XZY":
				this._z = Math.asin(-n(s, -1, 1)), Math.abs(s) < .9999999 ? (this._x = Math.atan2(d, c), this._y = Math.atan2(o, a)) : (this._x = Math.atan2(-h, p), this._y = 0);
				break;
			default:
				console.warn("THREE.Euler: .setFromRotationMatrix() encountered an unknown order: " + t)
		}
		return this._order = t, !1 !== i && this._onChangeCallback(), this
	}
	setFromQuaternion(e, t, i) {
		return _matrix.makeRotationFromQuaternion(e), this.setFromRotationMatrix(_matrix, t, i)
	}
	setFromVector3(e, t) {
		return this.set(e.x, e.y, e.z, t || this._order)
	}
	reorder(e) {
		return _quaternion$1.setFromEuler(this), this.setFromQuaternion(_quaternion$1, e)
	}
	equals(e) {
		return e._x === this._x && e._y === this._y && e._z === this._z && e._order === this._order
	}
	fromArray(e) {
		return this._x = e[0], this._y = e[1], this._z = e[2], void 0 !== e[3] && (this._order = e[3]), this._onChangeCallback(), this
	}
	toArray(e = [], t = 0) {
		return e[t] = this._x, e[t + 1] = this._y, e[t + 2] = this._z, e[t + 3] = this._order, e
	}
	toVector3(e) {
		return e ? e.set(this._x, this._y, this._z) : new Vector3(this._x, this._y, this._z)
	}
	_onChange(e) {
		return this._onChangeCallback = e, this
	}
	_onChangeCallback() {}
}
Euler.prototype.isEuler = !0, Euler.DefaultOrder = "XYZ", Euler.RotationOrders = ["XYZ", "YZX", "ZXY", "XZY", "YXZ", "ZYX"];
class Layers {
	constructor() {
		this.mask = 1
	}
	set(e) {
		this.mask = 1 << e | 0
	}
	enable(e) {
		this.mask |= 1 << e | 0
	}
	enableAll() {
		this.mask = -1
	}
	toggle(e) {
		this.mask ^= 1 << e | 0
	}
	disable(e) {
		this.mask &= ~(1 << e | 0)
	}
	disableAll() {
		this.mask = 0
	}
	test(e) {
		return 0 != (this.mask & e.mask)
	}
}
let _object3DId = 0;
const _v1$2 = new Vector3,
	_q1 = new Quaternion,
	_m1$1 = new Matrix4,
	_target = new Vector3,
	_position = new Vector3,
	_scale = new Vector3,
	_quaternion$2 = new Quaternion,
	_xAxis = new Vector3(1, 0, 0),
	_yAxis = new Vector3(0, 1, 0),
	_zAxis = new Vector3(0, 0, 1),
	_addedEvent = {
		type: "added"
	},
	_removedEvent = {
		type: "removed"
	};

function Object3D() {
	Object.defineProperty(this, "id", {
		value: _object3DId++
	}), this.uuid = MathUtils.generateUUID(), this.name = "", this.type = "Object3D", this.parent = null, this.children = [], this.up = Object3D.DefaultUp.clone();
	const e = new Vector3,
		t = new Euler,
		i = new Quaternion,
		n = new Vector3(1, 1, 1);
	t._onChange((function() {
		i.setFromEuler(t, !1)
	})), i._onChange((function() {
		t.setFromQuaternion(i, void 0, !1)
	})), Object.defineProperties(this, {
		position: {
			configurable: !0,
			enumerable: !0,
			value: e
		},
		rotation: {
			configurable: !0,
			enumerable: !0,
			value: t
		},
		quaternion: {
			configurable: !0,
			enumerable: !0,
			value: i
		},
		scale: {
			configurable: !0,
			enumerable: !0,
			value: n
		},
		modelViewMatrix: {
			value: new Matrix4
		},
		normalMatrix: {
			value: new Matrix3
		}
	}), this.matrix = new Matrix4, this.matrixWorld = new Matrix4, this.matrixAutoUpdate = Object3D.DefaultMatrixAutoUpdate, this.matrixWorldNeedsUpdate = !1, this.layers = new Layers, this.visible = !0, this.castShadow = !1, this.receiveShadow = !1, this.frustumCulled = !0, this.renderOrder = 0, this.animations = [], this.userData = {}
}
Object3D.DefaultUp = new Vector3(0, 1, 0), Object3D.DefaultMatrixAutoUpdate = !0, Object3D.prototype = Object.assign(Object.create(EventDispatcher.prototype), {
	constructor: Object3D,
	isObject3D: !0,
	onBeforeRender: function() {},
	onAfterRender: function() {},
	applyMatrix4: function(e) {
		this.matrixAutoUpdate && this.updateMatrix(), this.matrix.premultiply(e), this.matrix.decompose(this.position, this.quaternion, this.scale)
	},
	applyQuaternion: function(e) {
		return this.quaternion.premultiply(e), this
	},
	setRotationFromAxisAngle: function(e, t) {
		this.quaternion.setFromAxisAngle(e, t)
	},
	setRotationFromEuler: function(e) {
		this.quaternion.setFromEuler(e, !0)
	},
	setRotationFromMatrix: function(e) {
		this.quaternion.setFromRotationMatrix(e)
	},
	setRotationFromQuaternion: function(e) {
		this.quaternion.copy(e)
	},
	rotateOnAxis: function(e, t) {
		return _q1.setFromAxisAngle(e, t), this.quaternion.multiply(_q1), this
	},
	rotateOnWorldAxis: function(e, t) {
		return _q1.setFromAxisAngle(e, t), this.quaternion.premultiply(_q1), this
	},
	rotateX: function(e) {
		return this.rotateOnAxis(_xAxis, e)
	},
	rotateY: function(e) {
		return this.rotateOnAxis(_yAxis, e)
	},
	rotateZ: function(e) {
		return this.rotateOnAxis(_zAxis, e)
	},
	translateOnAxis: function(e, t) {
		return _v1$2.copy(e).applyQuaternion(this.quaternion), this.position.add(_v1$2.multiplyScalar(t)), this
	},
	translateX: function(e) {
		return this.translateOnAxis(_xAxis, e)
	},
	translateY: function(e) {
		return this.translateOnAxis(_yAxis, e)
	},
	translateZ: function(e) {
		return this.translateOnAxis(_zAxis, e)
	},
	localToWorld: function(e) {
		return e.applyMatrix4(this.matrixWorld)
	},
	worldToLocal: function(e) {
		return e.applyMatrix4(_m1$1.copy(this.matrixWorld).invert())
	},
	lookAt: function(e, t, i) {
		e.isVector3 ? _target.copy(e) : _target.set(e, t, i);
		const n = this.parent;
		this.updateWorldMatrix(!0, !1), _position.setFromMatrixPosition(this.matrixWorld), this.isCamera || this.isLight ? _m1$1.lookAt(_position, _target, this.up) : _m1$1.lookAt(_target, _position, this.up), this.quaternion.setFromRotationMatrix(_m1$1), n && (_m1$1.extractRotation(n.matrixWorld), _q1.setFromRotationMatrix(_m1$1), this.quaternion.premultiply(_q1.invert()))
	},
	add: function(e) {
		if (arguments.length > 1) {
			for (let e = 0; e < arguments.length; e++) this.add(arguments[e]);
			return this
		}
		return e === this ? (console.error("THREE.Object3D.add: object can't be added as a child of itself.", e), this) : (e && e.isObject3D ? (null !== e.parent && e.parent.remove(e), e.parent = this, this.children.push(e), e.dispatchEvent(_addedEvent)) : console.error("THREE.Object3D.add: object not an instance of THREE.Object3D.", e), this)
	},
	remove: function(e) {
		if (arguments.length > 1) {
			for (let e = 0; e < arguments.length; e++) this.remove(arguments[e]);
			return this
		}
		const t = this.children.indexOf(e);
		return -1 !== t && (e.parent = null, this.children.splice(t, 1), e.dispatchEvent(_removedEvent)), this
	},
	clear: function() {
		for (let e = 0; e < this.children.length; e++) {
			const t = this.children[e];
			t.parent = null, t.dispatchEvent(_removedEvent)
		}
		return this.children.length = 0, this
	},
	attach: function(e) {
		return this.updateWorldMatrix(!0, !1), _m1$1.copy(this.matrixWorld).invert(), null !== e.parent && (e.parent.updateWorldMatrix(!0, !1), _m1$1.multiply(e.parent.matrixWorld)), e.applyMatrix4(_m1$1), this.add(e), e.updateWorldMatrix(!1, !0), this
	},
	getObjectById: function(e) {
		return this.getObjectByProperty("id", e)
	},
	getObjectByName: function(e) {
		return this.getObjectByProperty("name", e)
	},
	getObjectByProperty: function(e, t) {
		if (this[e] === t) return this;
		for (let i = 0, n = this.children.length; i < n; i++) {
			const n = this.children[i].getObjectByProperty(e, t);
			if (void 0 !== n) return n
		}
	},
	getWorldPosition: function(e) {
		return void 0 === e && (console.warn("THREE.Object3D: .getWorldPosition() target is now required"), e = new Vector3), this.updateWorldMatrix(!0, !1), e.setFromMatrixPosition(this.matrixWorld)
	},
	getWorldQuaternion: function(e) {
		return void 0 === e && (console.warn("THREE.Object3D: .getWorldQuaternion() target is now required"), e = new Quaternion), this.updateWorldMatrix(!0, !1), this.matrixWorld.decompose(_position, e, _scale), e
	},
	getWorldScale: function(e) {
		return void 0 === e && (console.warn("THREE.Object3D: .getWorldScale() target is now required"), e = new Vector3), this.updateWorldMatrix(!0, !1), this.matrixWorld.decompose(_position, _quaternion$2, e), e
	},
	getWorldDirection: function(e) {
		void 0 === e && (console.warn("THREE.Object3D: .getWorldDirection() target is now required"), e = new Vector3), this.updateWorldMatrix(!0, !1);
		const t = this.matrixWorld.elements;
		return e.set(t[8], t[9], t[10]).normalize()
	},
	raycast: function() {},
	traverse: function(e) {
		e(this);
		const t = this.children;
		for (let i = 0, n = t.length; i < n; i++) t[i].traverse(e)
	},
	traverseVisible: function(e) {
		if (!1 === this.visible) return;
		e(this);
		const t = this.children;
		for (let i = 0, n = t.length; i < n; i++) t[i].traverseVisible(e)
	},
	traverseAncestors: function(e) {
		const t = this.parent;
		null !== t && (e(t), t.traverseAncestors(e))
	},
	updateMatrix: function() {
		this.matrix.compose(this.position, this.quaternion, this.scale), this.matrixWorldNeedsUpdate = !0
	},
	updateMatrixWorld: function(e) {
		this.matrixAutoUpdate && this.updateMatrix(), (this.matrixWorldNeedsUpdate || e) && (null === this.parent ? this.matrixWorld.copy(this.matrix) : this.matrixWorld.multiplyMatrices(this.parent.matrixWorld, this.matrix), this.matrixWorldNeedsUpdate = !1, e = !0);
		const t = this.children;
		for (let i = 0, n = t.length; i < n; i++) t[i].updateMatrixWorld(e)
	},
	updateWorldMatrix: function(e, t) {
		const i = this.parent;
		if (!0 === e && null !== i && i.updateWorldMatrix(!0, !1), this.matrixAutoUpdate && this.updateMatrix(), null === this.parent ? this.matrixWorld.copy(this.matrix) : this.matrixWorld.multiplyMatrices(this.parent.matrixWorld, this.matrix), !0 === t) {
			const e = this.children;
			for (let t = 0, i = e.length; t < i; t++) e[t].updateWorldMatrix(!1, !0)
		}
	},
	toJSON: function(e) {
		const t = void 0 === e || "string" == typeof e,
			i = {};
		t && (e = {
			geometries: {},
			materials: {},
			textures: {},
			images: {},
			shapes: {},
			skeletons: {},
			animations: {}
		}, i.metadata = {
			version: 4.5,
			type: "Object",
			generator: "Object3D.toJSON"
		});
		const n = {};

		function r(t, i) {
			return void 0 === t[i.uuid] && (t[i.uuid] = i.toJSON(e)), i.uuid
		}
		if (n.uuid = this.uuid, n.type = this.type, "" !== this.name && (n.name = this.name), !0 === this.castShadow && (n.castShadow = !0), !0 === this.receiveShadow && (n.receiveShadow = !0), !1 === this.visible && (n.visible = !1), !1 === this.frustumCulled && (n.frustumCulled = !1), 0 !== this.renderOrder && (n.renderOrder = this.renderOrder), "{}" !== JSON.stringify(this.userData) && (n.userData = this.userData), n.layers = this.layers.mask, n.matrix = this.matrix.toArray(), !1 === this.matrixAutoUpdate && (n.matrixAutoUpdate = !1), this.isInstancedMesh && (n.type = "InstancedMesh", n.count = this.count, n.instanceMatrix = this.instanceMatrix.toJSON()), this.isMesh || this.isLine || this.isPoints) {
			n.geometry = r(e.geometries, this.geometry);
			const t = this.geometry.parameters;
			if (void 0 !== t && void 0 !== t.shapes) {
				const i = t.shapes;
				if (Array.isArray(i))
					for (let t = 0, n = i.length; t < n; t++) {
						const n = i[t];
						r(e.shapes, n)
					} else r(e.shapes, i)
			}
		}
		if (this.isSkinnedMesh && (n.bindMode = this.bindMode, n.bindMatrix = this.bindMatrix.toArray(), void 0 !== this.skeleton && (r(e.skeletons, this.skeleton), n.skeleton = this.skeleton.uuid)), void 0 !== this.material)
			if (Array.isArray(this.material)) {
				const t = [];
				for (let i = 0, n = this.material.length; i < n; i++) t.push(r(e.materials, this.material[i]));
				n.material = t
			} else n.material = r(e.materials, this.material);
		if (this.children.length > 0) {
			n.children = [];
			for (let t = 0; t < this.children.length; t++) n.children.push(this.children[t].toJSON(e).object)
		}
		if (this.animations.length > 0) {
			n.animations = [];
			for (let t = 0; t < this.animations.length; t++) {
				const i = this.animations[t];
				n.animations.push(r(e.animations, i))
			}
		}
		if (t) {
			const t = a(e.geometries),
				n = a(e.materials),
				r = a(e.textures),
				s = a(e.images),
				o = a(e.shapes),
				l = a(e.skeletons),
				c = a(e.animations);
			t.length > 0 && (i.geometries = t), n.length > 0 && (i.materials = n), r.length > 0 && (i.textures = r), s.length > 0 && (i.images = s), o.length > 0 && (i.shapes = o), l.length > 0 && (i.skeletons = l), c.length > 0 && (i.animations = c)
		}
		return i.object = n, i;

		function a(e) {
			const t = [];
			for (const i in e) {
				const n = e[i];
				delete n.metadata, t.push(n)
			}
			return t
		}
	},
	clone: function(e) {
		return (new this.constructor).copy(this, e)
	},
	copy: function(e, t = !0) {
		if (this.name = e.name, this.up.copy(e.up), this.position.copy(e.position), this.rotation.order = e.rotation.order, this.quaternion.copy(e.quaternion), this.scale.copy(e.scale), this.matrix.copy(e.matrix), this.matrixWorld.copy(e.matrixWorld), this.matrixAutoUpdate = e.matrixAutoUpdate, this.matrixWorldNeedsUpdate = e.matrixWorldNeedsUpdate, this.layers.mask = e.layers.mask, this.visible = e.visible, this.castShadow = e.castShadow, this.receiveShadow = e.receiveShadow, this.frustumCulled = e.frustumCulled, this.renderOrder = e.renderOrder, this.userData = JSON.parse(JSON.stringify(e.userData)), !0 === t)
			for (let t = 0; t < e.children.length; t++) {
				const i = e.children[t];
				this.add(i.clone())
			}
		return this
	}
});
const _vector1 = new Vector3,
	_vector2 = new Vector3,
	_normalMatrix = new Matrix3;
class Plane {
	constructor(e = new Vector3(1, 0, 0), t = 0) {
		this.normal = e, this.constant = t
	}
	set(e, t) {
		return this.normal.copy(e), this.constant = t, this
	}
	setComponents(e, t, i, n) {
		return this.normal.set(e, t, i), this.constant = n, this
	}
	setFromNormalAndCoplanarPoint(e, t) {
		return this.normal.copy(e), this.constant = -t.dot(this.normal), this
	}
	setFromCoplanarPoints(e, t, i) {
		const n = _vector1.subVectors(i, t).cross(_vector2.subVectors(e, t)).normalize();
		return this.setFromNormalAndCoplanarPoint(n, e), this
	}
	copy(e) {
		return this.normal.copy(e.normal), this.constant = e.constant, this
	}
	normalize() {
		const e = 1 / this.normal.length();
		return this.normal.multiplyScalar(e), this.constant *= e, this
	}
	negate() {
		return this.constant *= -1, this.normal.negate(), this
	}
	distanceToPoint(e) {
		return this.normal.dot(e) + this.constant
	}
	distanceToSphere(e) {
		return this.distanceToPoint(e.center) - e.radius
	}
	projectPoint(e, t) {
		return void 0 === t && (console.warn("THREE.Plane: .projectPoint() target is now required"), t = new Vector3), t.copy(this.normal).multiplyScalar(-this.distanceToPoint(e)).add(e)
	}
	intersectLine(e, t) {
		void 0 === t && (console.warn("THREE.Plane: .intersectLine() target is now required"), t = new Vector3);
		const i = e.delta(_vector1),
			n = this.normal.dot(i);
		if (0 === n) return 0 === this.distanceToPoint(e.start) ? t.copy(e.start) : void 0;
		const r = -(e.start.dot(this.normal) + this.constant) / n;
		return r < 0 || r > 1 ? void 0 : t.copy(i).multiplyScalar(r).add(e.start)
	}
	intersectsLine(e) {
		const t = this.distanceToPoint(e.start),
			i = this.distanceToPoint(e.end);
		return t < 0 && i > 0 || i < 0 && t > 0
	}
	intersectsBox(e) {
		return e.intersectsPlane(this)
	}
	intersectsSphere(e) {
		return e.intersectsPlane(this)
	}
	coplanarPoint(e) {
		return void 0 === e && (console.warn("THREE.Plane: .coplanarPoint() target is now required"), e = new Vector3), e.copy(this.normal).multiplyScalar(-this.constant)
	}
	applyMatrix4(e, t) {
		const i = t || _normalMatrix.getNormalMatrix(e),
			n = this.coplanarPoint(_vector1).applyMatrix4(e),
			r = this.normal.applyMatrix3(i).normalize();
		return this.constant = -n.dot(r), this
	}
	translate(e) {
		return this.constant -= e.dot(this.normal), this
	}
	equals(e) {
		return e.normal.equals(this.normal) && e.constant === this.constant
	}
	clone() {
		return (new this.constructor).copy(this)
	}
}
Plane.prototype.isPlane = !0;
const _v0$1 = new Vector3,
	_v1$3 = new Vector3,
	_v2$1 = new Vector3,
	_v3 = new Vector3,
	_vab = new Vector3,
	_vac = new Vector3,
	_vbc = new Vector3,
	_vap = new Vector3,
	_vbp = new Vector3,
	_vcp = new Vector3;
class Triangle {
	constructor(e = new Vector3, t = new Vector3, i = new Vector3) {
		this.a = e, this.b = t, this.c = i
	}
	static getNormal(e, t, i, n) {
		void 0 === n && (console.warn("THREE.Triangle: .getNormal() target is now required"), n = new Vector3), n.subVectors(i, t), _v0$1.subVectors(e, t), n.cross(_v0$1);
		const r = n.lengthSq();
		return r > 0 ? n.multiplyScalar(1 / Math.sqrt(r)) : n.set(0, 0, 0)
	}
	static getBarycoord(e, t, i, n, r) {
		_v0$1.subVectors(n, t), _v1$3.subVectors(i, t), _v2$1.subVectors(e, t);
		const a = _v0$1.dot(_v0$1),
			s = _v0$1.dot(_v1$3),
			o = _v0$1.dot(_v2$1),
			l = _v1$3.dot(_v1$3),
			c = _v1$3.dot(_v2$1),
			h = a * l - s * s;
		if (void 0 === r && (console.warn("THREE.Triangle: .getBarycoord() target is now required"), r = new Vector3), 0 === h) return r.set(-2, -1, -1);
		const u = 1 / h,
			d = (l * o - s * c) * u,
			p = (a * c - s * o) * u;
		return r.set(1 - d - p, p, d)
	}
	static containsPoint(e, t, i, n) {
		return this.getBarycoord(e, t, i, n, _v3), _v3.x >= 0 && _v3.y >= 0 && _v3.x + _v3.y <= 1
	}
	static getUV(e, t, i, n, r, a, s, o) {
		return this.getBarycoord(e, t, i, n, _v3), o.set(0, 0), o.addScaledVector(r, _v3.x), o.addScaledVector(a, _v3.y), o.addScaledVector(s, _v3.z), o
	}
	static isFrontFacing(e, t, i, n) {
		return _v0$1.subVectors(i, t), _v1$3.subVectors(e, t), _v0$1.cross(_v1$3).dot(n) < 0
	}
	set(e, t, i) {
		return this.a.copy(e), this.b.copy(t), this.c.copy(i), this
	}
	setFromPointsAndIndices(e, t, i, n) {
		return this.a.copy(e[t]), this.b.copy(e[i]), this.c.copy(e[n]), this
	}
	clone() {
		return (new this.constructor).copy(this)
	}
	copy(e) {
		return this.a.copy(e.a), this.b.copy(e.b), this.c.copy(e.c), this
	}
	getArea() {
		return _v0$1.subVectors(this.c, this.b), _v1$3.subVectors(this.a, this.b), .5 * _v0$1.cross(_v1$3).length()
	}
	getMidpoint(e) {
		return void 0 === e && (console.warn("THREE.Triangle: .getMidpoint() target is now required"), e = new Vector3), e.addVectors(this.a, this.b).add(this.c).multiplyScalar(1 / 3)
	}
	getNormal(e) {
		return Triangle.getNormal(this.a, this.b, this.c, e)
	}
	getPlane(e) {
		return void 0 === e && (console.warn("THREE.Triangle: .getPlane() target is now required"), e = new Plane), e.setFromCoplanarPoints(this.a, this.b, this.c)
	}
	getBarycoord(e, t) {
		return Triangle.getBarycoord(e, this.a, this.b, this.c, t)
	}
	getUV(e, t, i, n, r) {
		return Triangle.getUV(e, this.a, this.b, this.c, t, i, n, r)
	}
	containsPoint(e) {
		return Triangle.containsPoint(e, this.a, this.b, this.c)
	}
	isFrontFacing(e) {
		return Triangle.isFrontFacing(this.a, this.b, this.c, e)
	}
	intersectsBox(e) {
		return e.intersectsTriangle(this)
	}
	closestPointToPoint(e, t) {
		void 0 === t && (console.warn("THREE.Triangle: .closestPointToPoint() target is now required"), t = new Vector3);
		const i = this.a,
			n = this.b,
			r = this.c;
		let a, s;
		_vab.subVectors(n, i), _vac.subVectors(r, i), _vap.subVectors(e, i);
		const o = _vab.dot(_vap),
			l = _vac.dot(_vap);
		if (o <= 0 && l <= 0) return t.copy(i);
		_vbp.subVectors(e, n);
		const c = _vab.dot(_vbp),
			h = _vac.dot(_vbp);
		if (c >= 0 && h <= c) return t.copy(n);
		const u = o * h - c * l;
		if (u <= 0 && o >= 0 && c <= 0) return a = o / (o - c), t.copy(i).addScaledVector(_vab, a);
		_vcp.subVectors(e, r);
		const d = _vab.dot(_vcp),
			p = _vac.dot(_vcp);
		if (p >= 0 && d <= p) return t.copy(r);
		const m = d * l - o * p;
		if (m <= 0 && l >= 0 && p <= 0) return s = l / (l - p), t.copy(i).addScaledVector(_vac, s);
		const A = c * p - d * h;
		if (A <= 0 && h - c >= 0 && d - p >= 0) return _vbc.subVectors(r, n), s = (h - c) / (h - c + (d - p)), t.copy(n).addScaledVector(_vbc, s);
		const g = 1 / (A + m + u);
		return a = m * g, s = u * g, t.copy(i).addScaledVector(_vab, a).addScaledVector(_vac, s)
	}
	equals(e) {
		return e.a.equals(this.a) && e.b.equals(this.b) && e.c.equals(this.c)
	}
}
let materialId = 0;

function Material$1() {
	Object.defineProperty(this, "id", {
		value: materialId++
	}), this.uuid = MathUtils.generateUUID(), this.name = "", this.type = "Material", this.fog = !0, this.blending = 1, this.side = 0, this.vertexColors = !1, this.opacity = 1, this.transparent = !1, this.blendSrc = 204, this.blendDst = 205, this.blendEquation = 100, this.blendSrcAlpha = null, this.blendDstAlpha = null, this.blendEquationAlpha = null, this.depthFunc = 3, this.depthTest = !0, this.depthWrite = !0, this.stencilWriteMask = 255, this.stencilFunc = 519, this.stencilRef = 0, this.stencilFuncMask = 255, this.stencilFail = 7680, this.stencilZFail = 7680, this.stencilZPass = 7680, this.stencilWrite = !1, this.clippingPlanes = null, this.clipIntersection = !1, this.clipShadows = !1, this.shadowSide = null, this.colorWrite = !0, this.precision = null, this.polygonOffset = !1, this.polygonOffsetFactor = 0, this.polygonOffsetUnits = 0, this.dithering = !1, this.alphaTest = 0, this.premultipliedAlpha = !1, this.visible = !0, this.toneMapped = !0, this.userData = {}, this.version = 0
}
Material$1.prototype = Object.assign(Object.create(EventDispatcher.prototype), {
	constructor: Material$1,
	isMaterial: !0,
	onBeforeCompile: function() {},
	customProgramCacheKey: function() {
		return this.onBeforeCompile.toString()
	},
	setValues: function(e) {
		if (void 0 !== e)
			for (const t in e) {
				const i = e[t];
				if (void 0 === i) {
					console.warn("THREE.Material: '" + t + "' parameter is undefined.");
					continue
				}
				if ("shading" === t) {
					console.warn("THREE." + this.type + ": .shading has been removed. Use the boolean .flatShading instead."), this.flatShading = 1 === i;
					continue
				}
				const n = this[t];
				void 0 !== n ? n && n.isColor ? n.set(i) : n && n.isVector3 && i && i.isVector3 ? n.copy(i) : this[t] = i : console.warn("THREE." + this.type + ": '" + t + "' is not a property of this material.")
			}
	},
	toJSON: function(e) {
		const t = void 0 === e || "string" == typeof e;
		t && (e = {
			textures: {},
			images: {}
		});
		const i = {
			metadata: {
				version: 4.5,
				type: "Material",
				generator: "Material.toJSON"
			}
		};

		function n(e) {
			const t = [];
			for (const i in e) {
				const n = e[i];
				delete n.metadata, t.push(n)
			}
			return t
		}
		if (i.uuid = this.uuid, i.type = this.type, "" !== this.name && (i.name = this.name), this.color && this.color.isColor && (i.color = this.color.getHex()), void 0 !== this.roughness && (i.roughness = this.roughness), void 0 !== this.metalness && (i.metalness = this.metalness), this.sheen && this.sheen.isColor && (i.sheen = this.sheen.getHex()), this.emissive && this.emissive.isColor && (i.emissive = this.emissive.getHex()), this.emissiveIntensity && 1 !== this.emissiveIntensity && (i.emissiveIntensity = this.emissiveIntensity), this.specular && this.specular.isColor && (i.specular = this.specular.getHex()), void 0 !== this.shininess && (i.shininess = this.shininess), void 0 !== this.clearcoat && (i.clearcoat = this.clearcoat), void 0 !== this.clearcoatRoughness && (i.clearcoatRoughness = this.clearcoatRoughness), this.clearcoatMap && this.clearcoatMap.isTexture && (i.clearcoatMap = this.clearcoatMap.toJSON(e).uuid), this.clearcoatRoughnessMap && this.clearcoatRoughnessMap.isTexture && (i.clearcoatRoughnessMap = this.clearcoatRoughnessMap.toJSON(e).uuid), this.clearcoatNormalMap && this.clearcoatNormalMap.isTexture && (i.clearcoatNormalMap = this.clearcoatNormalMap.toJSON(e).uuid, i.clearcoatNormalScale = this.clearcoatNormalScale.toArray()), this.map && this.map.isTexture && (i.map = this.map.toJSON(e).uuid), this.matcap && this.matcap.isTexture && (i.matcap = this.matcap.toJSON(e).uuid), this.alphaMap && this.alphaMap.isTexture && (i.alphaMap = this.alphaMap.toJSON(e).uuid), this.lightMap && this.lightMap.isTexture && (i.lightMap = this.lightMap.toJSON(e).uuid, i.lightMapIntensity = this.lightMapIntensity), this.aoMap && this.aoMap.isTexture && (i.aoMap = this.aoMap.toJSON(e).uuid, i.aoMapIntensity = this.aoMapIntensity), this.bumpMap && this.bumpMap.isTexture && (i.bumpMap = this.bumpMap.toJSON(e).uuid, i.bumpScale = this.bumpScale), this.normalMap && this.normalMap.isTexture && (i.normalMap = this.normalMap.toJSON(e).uuid, i.normalMapType = this.normalMapType, i.normalScale = this.normalScale.toArray()), this.displacementMap && this.displacementMap.isTexture && (i.displacementMap = this.displacementMap.toJSON(e).uuid, i.displacementScale = this.displacementScale, i.displacementBias = this.displacementBias), this.roughnessMap && this.roughnessMap.isTexture && (i.roughnessMap = this.roughnessMap.toJSON(e).uuid), this.metalnessMap && this.metalnessMap.isTexture && (i.metalnessMap = this.metalnessMap.toJSON(e).uuid), this.emissiveMap && this.emissiveMap.isTexture && (i.emissiveMap = this.emissiveMap.toJSON(e).uuid), this.specularMap && this.specularMap.isTexture && (i.specularMap = this.specularMap.toJSON(e).uuid), this.envMap && this.envMap.isTexture && (i.envMap = this.envMap.toJSON(e).uuid, i.reflectivity = this.reflectivity, i.refractionRatio = this.refractionRatio, void 0 !== this.combine && (i.combine = this.combine), void 0 !== this.envMapIntensity && (i.envMapIntensity = this.envMapIntensity)), this.gradientMap && this.gradientMap.isTexture && (i.gradientMap = this.gradientMap.toJSON(e).uuid), void 0 !== this.size && (i.size = this.size), void 0 !== this.sizeAttenuation && (i.sizeAttenuation = this.sizeAttenuation), 1 !== this.blending && (i.blending = this.blending), 0 !== this.side && (i.side = this.side), this.vertexColors && (i.vertexColors = !0), this.opacity < 1 && (i.opacity = this.opacity), !0 === this.transparent && (i.transparent = this.transparent), i.depthFunc = this.depthFunc, i.depthTest = this.depthTest, i.depthWrite = this.depthWrite, i.stencilWrite = this.stencilWrite, i.stencilWriteMask = this.stencilWriteMask, i.stencilFunc = this.stencilFunc, i.stencilRef = this.stencilRef, i.stencilFuncMask = this.stencilFuncMask, i.stencilFail = this.stencilFail, i.stencilZFail = this.stencilZFail, i.stencilZPass = this.stencilZPass, this.rotation && 0 !== this.rotation && (i.rotation = this.rotation), !0 === this.polygonOffset && (i.polygonOffset = !0), 0 !== this.polygonOffsetFactor && (i.polygonOffsetFactor = this.polygonOffsetFactor), 0 !== this.polygonOffsetUnits && (i.polygonOffsetUnits = this.polygonOffsetUnits), this.linewidth && 1 !== this.linewidth && (i.linewidth = this.linewidth), void 0 !== this.dashSize && (i.dashSize = this.dashSize), void 0 !== this.gapSize && (i.gapSize = this.gapSize), void 0 !== this.scale && (i.scale = this.scale), !0 === this.dithering && (i.dithering = !0), this.alphaTest > 0 && (i.alphaTest = this.alphaTest), !0 === this.premultipliedAlpha && (i.premultipliedAlpha = this.premultipliedAlpha), !0 === this.wireframe && (i.wireframe = this.wireframe), this.wireframeLinewidth > 1 && (i.wireframeLinewidth = this.wireframeLinewidth), "round" !== this.wireframeLinecap && (i.wireframeLinecap = this.wireframeLinecap), "round" !== this.wireframeLinejoin && (i.wireframeLinejoin = this.wireframeLinejoin), !0 === this.morphTargets && (i.morphTargets = !0), !0 === this.morphNormals && (i.morphNormals = !0), !0 === this.skinning && (i.skinning = !0), !0 === this.flatShading && (i.flatShading = this.flatShading), !1 === this.visible && (i.visible = !1), !1 === this.toneMapped && (i.toneMapped = !1), "{}" !== JSON.stringify(this.userData) && (i.userData = this.userData), t) {
			const t = n(e.textures),
				r = n(e.images);
			t.length > 0 && (i.textures = t), r.length > 0 && (i.images = r)
		}
		return i
	},
	clone: function() {
		return (new this.constructor).copy(this)
	},
	copy: function(e) {
		this.name = e.name, this.fog = e.fog, this.blending = e.blending, this.side = e.side, this.vertexColors = e.vertexColors, this.opacity = e.opacity, this.transparent = e.transparent, this.blendSrc = e.blendSrc, this.blendDst = e.blendDst, this.blendEquation = e.blendEquation, this.blendSrcAlpha = e.blendSrcAlpha, this.blendDstAlpha = e.blendDstAlpha, this.blendEquationAlpha = e.blendEquationAlpha, this.depthFunc = e.depthFunc, this.depthTest = e.depthTest, this.depthWrite = e.depthWrite, this.stencilWriteMask = e.stencilWriteMask, this.stencilFunc = e.stencilFunc, this.stencilRef = e.stencilRef, this.stencilFuncMask = e.stencilFuncMask, this.stencilFail = e.stencilFail, this.stencilZFail = e.stencilZFail, this.stencilZPass = e.stencilZPass, this.stencilWrite = e.stencilWrite;
		const t = e.clippingPlanes;
		let i = null;
		if (null !== t) {
			const e = t.length;
			i = new Array(e);
			for (let n = 0; n !== e; ++n) i[n] = t[n].clone()
		}
		return this.clippingPlanes = i, this.clipIntersection = e.clipIntersection, this.clipShadows = e.clipShadows, this.shadowSide = e.shadowSide, this.colorWrite = e.colorWrite, this.precision = e.precision, this.polygonOffset = e.polygonOffset, this.polygonOffsetFactor = e.polygonOffsetFactor, this.polygonOffsetUnits = e.polygonOffsetUnits, this.dithering = e.dithering, this.alphaTest = e.alphaTest, this.premultipliedAlpha = e.premultipliedAlpha, this.visible = e.visible, this.toneMapped = e.toneMapped, this.userData = JSON.parse(JSON.stringify(e.userData)), this
	},
	dispose: function() {
		this.dispatchEvent({
			type: "dispose"
		})
	}
}), Object.defineProperty(Material$1.prototype, "needsUpdate", {
	set: function(e) {
		!0 === e && this.version++
	}
});
const _colorKeywords = {
		aliceblue: 15792383,
		antiquewhite: 16444375,
		aqua: 65535,
		aquamarine: 8388564,
		azure: 15794175,
		beige: 16119260,
		bisque: 16770244,
		black: 0,
		blanchedalmond: 16772045,
		blue: 255,
		blueviolet: 9055202,
		brown: 10824234,
		burlywood: 14596231,
		cadetblue: 6266528,
		chartreuse: 8388352,
		chocolate: 13789470,
		coral: 16744272,
		cornflowerblue: 6591981,
		cornsilk: 16775388,
		crimson: 14423100,
		cyan: 65535,
		darkblue: 139,
		darkcyan: 35723,
		darkgoldenrod: 12092939,
		darkgray: 11119017,
		darkgreen: 25600,
		darkgrey: 11119017,
		darkkhaki: 12433259,
		darkmagenta: 9109643,
		darkolivegreen: 5597999,
		darkorange: 16747520,
		darkorchid: 10040012,
		darkred: 9109504,
		darksalmon: 15308410,
		darkseagreen: 9419919,
		darkslateblue: 4734347,
		darkslategray: 3100495,
		darkslategrey: 3100495,
		darkturquoise: 52945,
		darkviolet: 9699539,
		deeppink: 16716947,
		deepskyblue: 49151,
		dimgray: 6908265,
		dimgrey: 6908265,
		dodgerblue: 2003199,
		firebrick: 11674146,
		floralwhite: 16775920,
		forestgreen: 2263842,
		fuchsia: 16711935,
		gainsboro: 14474460,
		ghostwhite: 16316671,
		gold: 16766720,
		goldenrod: 14329120,
		gray: 8421504,
		green: 32768,
		greenyellow: 11403055,
		grey: 8421504,
		honeydew: 15794160,
		hotpink: 16738740,
		indianred: 13458524,
		indigo: 4915330,
		ivory: 16777200,
		khaki: 15787660,
		lavender: 15132410,
		lavenderblush: 16773365,
		lawngreen: 8190976,
		lemonchiffon: 16775885,
		lightblue: 11393254,
		lightcoral: 15761536,
		lightcyan: 14745599,
		lightgoldenrodyellow: 16448210,
		lightgray: 13882323,
		lightgreen: 9498256,
		lightgrey: 13882323,
		lightpink: 16758465,
		lightsalmon: 16752762,
		lightseagreen: 2142890,
		lightskyblue: 8900346,
		lightslategray: 7833753,
		lightslategrey: 7833753,
		lightsteelblue: 11584734,
		lightyellow: 16777184,
		lime: 65280,
		limegreen: 3329330,
		linen: 16445670,
		magenta: 16711935,
		maroon: 8388608,
		mediumaquamarine: 6737322,
		mediumblue: 205,
		mediumorchid: 12211667,
		mediumpurple: 9662683,
		mediumseagreen: 3978097,
		mediumslateblue: 8087790,
		mediumspringgreen: 64154,
		mediumturquoise: 4772300,
		mediumvioletred: 13047173,
		midnightblue: 1644912,
		mintcream: 16121850,
		mistyrose: 16770273,
		moccasin: 16770229,
		navajowhite: 16768685,
		navy: 128,
		oldlace: 16643558,
		olive: 8421376,
		olivedrab: 7048739,
		orange: 16753920,
		orangered: 16729344,
		orchid: 14315734,
		palegoldenrod: 15657130,
		palegreen: 10025880,
		paleturquoise: 11529966,
		palevioletred: 14381203,
		papayawhip: 16773077,
		peachpuff: 16767673,
		peru: 13468991,
		pink: 16761035,
		plum: 14524637,
		powderblue: 11591910,
		purple: 8388736,
		rebeccapurple: 6697881,
		red: 16711680,
		rosybrown: 12357519,
		royalblue: 4286945,
		saddlebrown: 9127187,
		salmon: 16416882,
		sandybrown: 16032864,
		seagreen: 3050327,
		seashell: 16774638,
		sienna: 10506797,
		silver: 12632256,
		skyblue: 8900331,
		slateblue: 6970061,
		slategray: 7372944,
		slategrey: 7372944,
		snow: 16775930,
		springgreen: 65407,
		steelblue: 4620980,
		tan: 13808780,
		teal: 32896,
		thistle: 14204888,
		tomato: 16737095,
		turquoise: 4251856,
		violet: 15631086,
		wheat: 16113331,
		white: 16777215,
		whitesmoke: 16119285,
		yellow: 16776960,
		yellowgreen: 10145074
	},
	_hslA = {
		h: 0,
		s: 0,
		l: 0
	},
	_hslB = {
		h: 0,
		s: 0,
		l: 0
	};

function hue2rgb(e, t, i) {
	return i < 0 && (i += 1), i > 1 && (i -= 1), i < 1 / 6 ? e + 6 * (t - e) * i : i < .5 ? t : i < 2 / 3 ? e + 6 * (t - e) * (2 / 3 - i) : e
}

function SRGBToLinear(e) {
	return e < .04045 ? .0773993808 * e : Math.pow(.9478672986 * e + .0521327014, 2.4)
}

function LinearToSRGB(e) {
	return e < .0031308 ? 12.92 * e : 1.055 * Math.pow(e, .41666) - .055
}
class Color {
	constructor(e, t, i) {
		return void 0 === t && void 0 === i ? this.set(e) : this.setRGB(e, t, i)
	}
	set(e) {
		return e && e.isColor ? this.copy(e) : "number" == typeof e ? this.setHex(e) : "string" == typeof e && this.setStyle(e), this
	}
	setScalar(e) {
		return this.r = e, this.g = e, this.b = e, this
	}
	setHex(e) {
		return e = Math.floor(e), this.r = (e >> 16 & 255) / 255, this.g = (e >> 8 & 255) / 255, this.b = (255 & e) / 255, this
	}
	setRGB(e, t, i) {
		return this.r = e, this.g = t, this.b = i, this
	}
	setHSL(e, t, i) {
		if (e = MathUtils.euclideanModulo(e, 1), t = MathUtils.clamp(t, 0, 1), i = MathUtils.clamp(i, 0, 1), 0 === t) this.r = this.g = this.b = i;
		else {
			const n = i <= .5 ? i * (1 + t) : i + t - i * t,
				r = 2 * i - n;
			this.r = hue2rgb(r, n, e + 1 / 3), this.g = hue2rgb(r, n, e), this.b = hue2rgb(r, n, e - 1 / 3)
		}
		return this
	}
	setStyle(e) {
		function t(t) {
			void 0 !== t && parseFloat(t) < 1 && console.warn("THREE.Color: Alpha component of " + e + " will be ignored.")
		}
		let i;
		if (i = /^((?:rgb|hsl)a?)\(([^\)]*)\)/.exec(e)) {
			let e;
			const n = i[1],
				r = i[2];
			switch (n) {
				case "rgb":
				case "rgba":
					if (e = /^\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*(?:,\s*(\d*\.?\d+)\s*)?$/.exec(r)) return this.r = Math.min(255, parseInt(e[1], 10)) / 255, this.g = Math.min(255, parseInt(e[2], 10)) / 255, this.b = Math.min(255, parseInt(e[3], 10)) / 255, t(e[4]), this;
					if (e = /^\s*(\d+)\%\s*,\s*(\d+)\%\s*,\s*(\d+)\%\s*(?:,\s*(\d*\.?\d+)\s*)?$/.exec(r)) return this.r = Math.min(100, parseInt(e[1], 10)) / 100, this.g = Math.min(100, parseInt(e[2], 10)) / 100, this.b = Math.min(100, parseInt(e[3], 10)) / 100, t(e[4]), this;
					break;
				case "hsl":
				case "hsla":
					if (e = /^\s*(\d*\.?\d+)\s*,\s*(\d+)\%\s*,\s*(\d+)\%\s*(?:,\s*(\d*\.?\d+)\s*)?$/.exec(r)) {
						const i = parseFloat(e[1]) / 360,
							n = parseInt(e[2], 10) / 100,
							r = parseInt(e[3], 10) / 100;
						return t(e[4]), this.setHSL(i, n, r)
					}
			}
		} else if (i = /^\#([A-Fa-f\d]+)$/.exec(e)) {
			const e = i[1],
				t = e.length;
			if (3 === t) return this.r = parseInt(e.charAt(0) + e.charAt(0), 16) / 255, this.g = parseInt(e.charAt(1) + e.charAt(1), 16) / 255, this.b = parseInt(e.charAt(2) + e.charAt(2), 16) / 255, this;
			if (6 === t) return this.r = parseInt(e.charAt(0) + e.charAt(1), 16) / 255, this.g = parseInt(e.charAt(2) + e.charAt(3), 16) / 255, this.b = parseInt(e.charAt(4) + e.charAt(5), 16) / 255, this
		}
		return e && e.length > 0 ? this.setColorName(e) : this
	}
	setColorName(e) {
		const t = _colorKeywords[e];
		return void 0 !== t ? this.setHex(t) : console.warn("THREE.Color: Unknown color " + e), this
	}
	clone() {
		return new this.constructor(this.r, this.g, this.b)
	}
	copy(e) {
		return this.r = e.r, this.g = e.g, this.b = e.b, this
	}
	copyGammaToLinear(e, t = 2) {
		return this.r = Math.pow(e.r, t), this.g = Math.pow(e.g, t), this.b = Math.pow(e.b, t), this
	}
	copyLinearToGamma(e, t = 2) {
		const i = t > 0 ? 1 / t : 1;
		return this.r = Math.pow(e.r, i), this.g = Math.pow(e.g, i), this.b = Math.pow(e.b, i), this
	}
	convertGammaToLinear(e) {
		return this.copyGammaToLinear(this, e), this
	}
	convertLinearToGamma(e) {
		return this.copyLinearToGamma(this, e), this
	}
	copySRGBToLinear(e) {
		return this.r = SRGBToLinear(e.r), this.g = SRGBToLinear(e.g), this.b = SRGBToLinear(e.b), this
	}
	copyLinearToSRGB(e) {
		return this.r = LinearToSRGB(e.r), this.g = LinearToSRGB(e.g), this.b = LinearToSRGB(e.b), this
	}
	convertSRGBToLinear() {
		return this.copySRGBToLinear(this), this
	}
	convertLinearToSRGB() {
		return this.copyLinearToSRGB(this), this
	}
	getHex() {
		return 255 * this.r << 16 ^ 255 * this.g << 8 ^ 255 * this.b << 0
	}
	getHexString() {
		return ("000000" + this.getHex().toString(16)).slice(-6)
	}
	getHSL(e) {
		void 0 === e && (console.warn("THREE.Color: .getHSL() target is now required"), e = {
			h: 0,
			s: 0,
			l: 0
		});
		const t = this.r,
			i = this.g,
			n = this.b,
			r = Math.max(t, i, n),
			a = Math.min(t, i, n);
		let s, o;
		const l = (a + r) / 2;
		if (a === r) s = 0, o = 0;
		else {
			const e = r - a;
			switch (o = l <= .5 ? e / (r + a) : e / (2 - r - a), r) {
				case t:
					s = (i - n) / e + (i < n ? 6 : 0);
					break;
				case i:
					s = (n - t) / e + 2;
					break;
				case n:
					s = (t - i) / e + 4
			}
			s /= 6
		}
		return e.h = s, e.s = o, e.l = l, e
	}
	getStyle() {
		return "rgb(" + (255 * this.r | 0) + "," + (255 * this.g | 0) + "," + (255 * this.b | 0) + ")"
	}
	offsetHSL(e, t, i) {
		return this.getHSL(_hslA), _hslA.h += e, _hslA.s += t, _hslA.l += i, this.setHSL(_hslA.h, _hslA.s, _hslA.l), this
	}
	add(e) {
		return this.r += e.r, this.g += e.g, this.b += e.b, this
	}
	addColors(e, t) {
		return this.r = e.r + t.r, this.g = e.g + t.g, this.b = e.b + t.b, this
	}
	addScalar(e) {
		return this.r += e, this.g += e, this.b += e, this
	}
	sub(e) {
		return this.r = Math.max(0, this.r - e.r), this.g = Math.max(0, this.g - e.g), this.b = Math.max(0, this.b - e.b), this
	}
	multiply(e) {
		return this.r *= e.r, this.g *= e.g, this.b *= e.b, this
	}
	multiplyScalar(e) {
		return this.r *= e, this.g *= e, this.b *= e, this
	}
	lerp(e, t) {
		return this.r += (e.r - this.r) * t, this.g += (e.g - this.g) * t, this.b += (e.b - this.b) * t, this
	}
	lerpColors(e, t, i) {
		return this.r = e.r + (t.r - e.r) * i, this.g = e.g + (t.g - e.g) * i, this.b = e.b + (t.b - e.b) * i, this
	}
	lerpHSL(e, t) {
		this.getHSL(_hslA), e.getHSL(_hslB);
		const i = MathUtils.lerp(_hslA.h, _hslB.h, t),
			n = MathUtils.lerp(_hslA.s, _hslB.s, t),
			r = MathUtils.lerp(_hslA.l, _hslB.l, t);
		return this.setHSL(i, n, r), this
	}
	equals(e) {
		return e.r === this.r && e.g === this.g && e.b === this.b
	}
	fromArray(e, t = 0) {
		return this.r = e[t], this.g = e[t + 1], this.b = e[t + 2], this
	}
	toArray(e = [], t = 0) {
		return e[t] = this.r, e[t + 1] = this.g, e[t + 2] = this.b, e
	}
	fromBufferAttribute(e, t) {
		return this.r = e.getX(t), this.g = e.getY(t), this.b = e.getZ(t), !0 === e.normalized && (this.r /= 255, this.g /= 255, this.b /= 255), this
	}
	toJSON() {
		return this.getHex()
	}
}
Color.NAMES = _colorKeywords, Color.prototype.isColor = !0, Color.prototype.r = 1, Color.prototype.g = 1, Color.prototype.b = 1;
class MeshBasicMaterial extends Material$1 {
	constructor(e) {
		super(), this.type = "MeshBasicMaterial", this.color = new Color(16777215), this.map = null, this.lightMap = null, this.lightMapIntensity = 1, this.aoMap = null, this.aoMapIntensity = 1, this.specularMap = null, this.alphaMap = null, this.envMap = null, this.combine = 0, this.reflectivity = 1, this.refractionRatio = .98, this.wireframe = !1, this.wireframeLinewidth = 1, this.wireframeLinecap = "round", this.wireframeLinejoin = "round", this.skinning = !1, this.morphTargets = !1, this.setValues(e)
	}
	copy(e) {
		return super.copy(e), this.color.copy(e.color), this.map = e.map, this.lightMap = e.lightMap, this.lightMapIntensity = e.lightMapIntensity, this.aoMap = e.aoMap, this.aoMapIntensity = e.aoMapIntensity, this.specularMap = e.specularMap, this.alphaMap = e.alphaMap, this.envMap = e.envMap, this.combine = e.combine, this.reflectivity = e.reflectivity, this.refractionRatio = e.refractionRatio, this.wireframe = e.wireframe, this.wireframeLinewidth = e.wireframeLinewidth, this.wireframeLinecap = e.wireframeLinecap, this.wireframeLinejoin = e.wireframeLinejoin, this.skinning = e.skinning, this.morphTargets = e.morphTargets, this
	}
}
MeshBasicMaterial.prototype.isMeshBasicMaterial = !0;
const _vector$3 = new Vector3,
	_vector2$1 = new Vector2;

function BufferAttribute(e, t, i) {
	if (Array.isArray(e)) throw new TypeError("THREE.BufferAttribute: array should be a Typed Array.");
	this.name = "", this.array = e, this.itemSize = t, this.count = void 0 !== e ? e.length / t : 0, this.normalized = !0 === i, this.usage = 35044, this.updateRange = {
		offset: 0,
		count: -1
	}, this.version = 0
}

function Int8BufferAttribute(e, t, i) {
	BufferAttribute.call(this, new Int8Array(e), t, i)
}

function Uint8BufferAttribute(e, t, i) {
	BufferAttribute.call(this, new Uint8Array(e), t, i)
}

function Uint8ClampedBufferAttribute(e, t, i) {
	BufferAttribute.call(this, new Uint8ClampedArray(e), t, i)
}

function Int16BufferAttribute(e, t, i) {
	BufferAttribute.call(this, new Int16Array(e), t, i)
}

function Uint16BufferAttribute(e, t, i) {
	BufferAttribute.call(this, new Uint16Array(e), t, i)
}

function Int32BufferAttribute(e, t, i) {
	BufferAttribute.call(this, new Int32Array(e), t, i)
}

function Uint32BufferAttribute(e, t, i) {
	BufferAttribute.call(this, new Uint32Array(e), t, i)
}

function Float16BufferAttribute(e, t, i) {
	BufferAttribute.call(this, new Uint16Array(e), t, i)
}

function Float32BufferAttribute(e, t, i) {
	BufferAttribute.call(this, new Float32Array(e), t, i)
}

function Float64BufferAttribute(e, t, i) {
	BufferAttribute.call(this, new Float64Array(e), t, i)
}

function arrayMax(e) {
	if (0 === e.length) return -1 / 0;
	let t = e[0];
	for (let i = 1, n = e.length; i < n; ++i) e[i] > t && (t = e[i]);
	return t
}
Object.defineProperty(BufferAttribute.prototype, "needsUpdate", {
	set: function(e) {
		!0 === e && this.version++
	}
}), Object.assign(BufferAttribute.prototype, {
	isBufferAttribute: !0,
	onUploadCallback: function() {},
	setUsage: function(e) {
		return this.usage = e, this
	},
	copy: function(e) {
		return this.name = e.name, this.array = new e.array.constructor(e.array), this.itemSize = e.itemSize, this.count = e.count, this.normalized = e.normalized, this.usage = e.usage, this
	},
	copyAt: function(e, t, i) {
		e *= this.itemSize, i *= t.itemSize;
		for (let n = 0, r = this.itemSize; n < r; n++) this.array[e + n] = t.array[i + n];
		return this
	},
	copyArray: function(e) {
		return this.array.set(e), this
	},
	copyColorsArray: function(e) {
		const t = this.array;
		let i = 0;
		for (let n = 0, r = e.length; n < r; n++) {
			let r = e[n];
			void 0 === r && (console.warn("THREE.BufferAttribute.copyColorsArray(): color is undefined", n), r = new Color), t[i++] = r.r, t[i++] = r.g, t[i++] = r.b
		}
		return this
	},
	copyVector2sArray: function(e) {
		const t = this.array;
		let i = 0;
		for (let n = 0, r = e.length; n < r; n++) {
			let r = e[n];
			void 0 === r && (console.warn("THREE.BufferAttribute.copyVector2sArray(): vector is undefined", n), r = new Vector2), t[i++] = r.x, t[i++] = r.y
		}
		return this
	},
	copyVector3sArray: function(e) {
		const t = this.array;
		let i = 0;
		for (let n = 0, r = e.length; n < r; n++) {
			let r = e[n];
			void 0 === r && (console.warn("THREE.BufferAttribute.copyVector3sArray(): vector is undefined", n), r = new Vector3), t[i++] = r.x, t[i++] = r.y, t[i++] = r.z
		}
		return this
	},
	copyVector4sArray: function(e) {
		const t = this.array;
		let i = 0;
		for (let n = 0, r = e.length; n < r; n++) {
			let r = e[n];
			void 0 === r && (console.warn("THREE.BufferAttribute.copyVector4sArray(): vector is undefined", n), r = new Vector4), t[i++] = r.x, t[i++] = r.y, t[i++] = r.z, t[i++] = r.w
		}
		return this
	},
	applyMatrix3: function(e) {
		if (2 === this.itemSize)
			for (let t = 0, i = this.count; t < i; t++) _vector2$1.fromBufferAttribute(this, t), _vector2$1.applyMatrix3(e), this.setXY(t, _vector2$1.x, _vector2$1.y);
		else if (3 === this.itemSize)
			for (let t = 0, i = this.count; t < i; t++) _vector$3.fromBufferAttribute(this, t), _vector$3.applyMatrix3(e), this.setXYZ(t, _vector$3.x, _vector$3.y, _vector$3.z);
		return this
	},
	applyMatrix4: function(e) {
		for (let t = 0, i = this.count; t < i; t++) _vector$3.x = this.getX(t), _vector$3.y = this.getY(t), _vector$3.z = this.getZ(t), _vector$3.applyMatrix4(e), this.setXYZ(t, _vector$3.x, _vector$3.y, _vector$3.z);
		return this
	},
	applyNormalMatrix: function(e) {
		for (let t = 0, i = this.count; t < i; t++) _vector$3.x = this.getX(t), _vector$3.y = this.getY(t), _vector$3.z = this.getZ(t), _vector$3.applyNormalMatrix(e), this.setXYZ(t, _vector$3.x, _vector$3.y, _vector$3.z);
		return this
	},
	transformDirection: function(e) {
		for (let t = 0, i = this.count; t < i; t++) _vector$3.x = this.getX(t), _vector$3.y = this.getY(t), _vector$3.z = this.getZ(t), _vector$3.transformDirection(e), this.setXYZ(t, _vector$3.x, _vector$3.y, _vector$3.z);
		return this
	},
	set: function(e, t = 0) {
		return this.array.set(e, t), this
	},
	getX: function(e) {
		return this.array[e * this.itemSize]
	},
	setX: function(e, t) {
		return this.array[e * this.itemSize] = t, this
	},
	getY: function(e) {
		return this.array[e * this.itemSize + 1]
	},
	setY: function(e, t) {
		return this.array[e * this.itemSize + 1] = t, this
	},
	getZ: function(e) {
		return this.array[e * this.itemSize + 2]
	},
	setZ: function(e, t) {
		return this.array[e * this.itemSize + 2] = t, this
	},
	getW: function(e) {
		return this.array[e * this.itemSize + 3]
	},
	setW: function(e, t) {
		return this.array[e * this.itemSize + 3] = t, this
	},
	setXY: function(e, t, i) {
		return e *= this.itemSize, this.array[e + 0] = t, this.array[e + 1] = i, this
	},
	setXYZ: function(e, t, i, n) {
		return e *= this.itemSize, this.array[e + 0] = t, this.array[e + 1] = i, this.array[e + 2] = n, this
	},
	setXYZW: function(e, t, i, n, r) {
		return e *= this.itemSize, this.array[e + 0] = t, this.array[e + 1] = i, this.array[e + 2] = n, this.array[e + 3] = r, this
	},
	onUpload: function(e) {
		return this.onUploadCallback = e, this
	},
	clone: function() {
		return new this.constructor(this.array, this.itemSize).copy(this)
	},
	toJSON: function() {
		return {
			itemSize: this.itemSize,
			type: this.array.constructor.name,
			array: Array.prototype.slice.call(this.array),
			normalized: this.normalized
		}
	}
}), Int8BufferAttribute.prototype = Object.create(BufferAttribute.prototype), Int8BufferAttribute.prototype.constructor = Int8BufferAttribute, Uint8BufferAttribute.prototype = Object.create(BufferAttribute.prototype), Uint8BufferAttribute.prototype.constructor = Uint8BufferAttribute, Uint8ClampedBufferAttribute.prototype = Object.create(BufferAttribute.prototype), Uint8ClampedBufferAttribute.prototype.constructor = Uint8ClampedBufferAttribute, Int16BufferAttribute.prototype = Object.create(BufferAttribute.prototype), Int16BufferAttribute.prototype.constructor = Int16BufferAttribute, Uint16BufferAttribute.prototype = Object.create(BufferAttribute.prototype), Uint16BufferAttribute.prototype.constructor = Uint16BufferAttribute, Int32BufferAttribute.prototype = Object.create(BufferAttribute.prototype), Int32BufferAttribute.prototype.constructor = Int32BufferAttribute, Uint32BufferAttribute.prototype = Object.create(BufferAttribute.prototype), Uint32BufferAttribute.prototype.constructor = Uint32BufferAttribute, Float16BufferAttribute.prototype = Object.create(BufferAttribute.prototype), Float16BufferAttribute.prototype.constructor = Float16BufferAttribute, Float16BufferAttribute.prototype.isFloat16BufferAttribute = !0, Float32BufferAttribute.prototype = Object.create(BufferAttribute.prototype), Float32BufferAttribute.prototype.constructor = Float32BufferAttribute, Float64BufferAttribute.prototype = Object.create(BufferAttribute.prototype), Float64BufferAttribute.prototype.constructor = Float64BufferAttribute;
let _id = 0;
const _m1$2 = new Matrix4,
	_obj = new Object3D,
	_offset = new Vector3,
	_box$2 = new Box3,
	_boxMorphTargets = new Box3,
	_vector$4 = new Vector3;

function BufferGeometry() {
	Object.defineProperty(this, "id", {
		value: _id++
	}), this.uuid = MathUtils.generateUUID(), this.name = "", this.type = "BufferGeometry", this.index = null, this.attributes = {}, this.morphAttributes = {}, this.morphTargetsRelative = !1, this.groups = [], this.boundingBox = null, this.boundingSphere = null, this.drawRange = {
		start: 0,
		count: 1 / 0
	}, this.userData = {}
}
BufferGeometry.prototype = Object.assign(Object.create(EventDispatcher.prototype), {
	constructor: BufferGeometry,
	isBufferGeometry: !0,
	getIndex: function() {
		return this.index
	},
	setIndex: function(e) {
		return Array.isArray(e) ? this.index = new(arrayMax(e) > 65535 ? Uint32BufferAttribute : Uint16BufferAttribute)(e, 1) : this.index = e, this
	},
	getAttribute: function(e) {
		return this.attributes[e]
	},
	setAttribute: function(e, t) {
		return this.attributes[e] = t, this
	},
	deleteAttribute: function(e) {
		return delete this.attributes[e], this
	},
	hasAttribute: function(e) {
		return void 0 !== this.attributes[e]
	},
	addGroup: function(e, t, i = 0) {
		this.groups.push({
			start: e,
			count: t,
			materialIndex: i
		})
	},
	clearGroups: function() {
		this.groups = []
	},
	setDrawRange: function(e, t) {
		this.drawRange.start = e, this.drawRange.count = t
	},
	applyMatrix4: function(e) {
		const t = this.attributes.position;
		void 0 !== t && (t.applyMatrix4(e), t.needsUpdate = !0);
		const i = this.attributes.normal;
		if (void 0 !== i) {
			const t = (new Matrix3).getNormalMatrix(e);
			i.applyNormalMatrix(t), i.needsUpdate = !0
		}
		const n = this.attributes.tangent;
		return void 0 !== n && (n.transformDirection(e), n.needsUpdate = !0), null !== this.boundingBox && this.computeBoundingBox(), null !== this.boundingSphere && this.computeBoundingSphere(), this
	},
	rotateX: function(e) {
		return _m1$2.makeRotationX(e), this.applyMatrix4(_m1$2), this
	},
	rotateY: function(e) {
		return _m1$2.makeRotationY(e), this.applyMatrix4(_m1$2), this
	},
	rotateZ: function(e) {
		return _m1$2.makeRotationZ(e), this.applyMatrix4(_m1$2), this
	},
	translate: function(e, t, i) {
		return _m1$2.makeTranslation(e, t, i), this.applyMatrix4(_m1$2), this
	},
	scale: function(e, t, i) {
		return _m1$2.makeScale(e, t, i), this.applyMatrix4(_m1$2), this
	},
	lookAt: function(e) {
		return _obj.lookAt(e), _obj.updateMatrix(), this.applyMatrix4(_obj.matrix), this
	},
	center: function() {
		return this.computeBoundingBox(), this.boundingBox.getCenter(_offset).negate(), this.translate(_offset.x, _offset.y, _offset.z), this
	},
	setFromPoints: function(e) {
		const t = [];
		for (let i = 0, n = e.length; i < n; i++) {
			const n = e[i];
			t.push(n.x, n.y, n.z || 0)
		}
		return this.setAttribute("position", new Float32BufferAttribute(t, 3)), this
	},
	computeBoundingBox: function() {
		null === this.boundingBox && (this.boundingBox = new Box3);
		const e = this.attributes.position,
			t = this.morphAttributes.position;
		if (e && e.isGLBufferAttribute) return console.error('THREE.BufferGeometry.computeBoundingBox(): GLBufferAttribute requires a manual bounding box. Alternatively set "mesh.frustumCulled" to "false".', this), void this.boundingBox.set(new Vector3(-1 / 0, -1 / 0, -1 / 0), new Vector3(1 / 0, 1 / 0, 1 / 0));
		if (void 0 !== e) {
			if (this.boundingBox.setFromBufferAttribute(e), t)
				for (let e = 0, i = t.length; e < i; e++) {
					const i = t[e];
					_box$2.setFromBufferAttribute(i), this.morphTargetsRelative ? (_vector$4.addVectors(this.boundingBox.min, _box$2.min), this.boundingBox.expandByPoint(_vector$4), _vector$4.addVectors(this.boundingBox.max, _box$2.max), this.boundingBox.expandByPoint(_vector$4)) : (this.boundingBox.expandByPoint(_box$2.min), this.boundingBox.expandByPoint(_box$2.max))
				}
		} else this.boundingBox.makeEmpty();
		(isNaN(this.boundingBox.min.x) || isNaN(this.boundingBox.min.y) || isNaN(this.boundingBox.min.z)) && console.error('THREE.BufferGeometry.computeBoundingBox(): Computed min/max have NaN values. The "position" attribute is likely to have NaN values.', this)
	},
	computeBoundingSphere: function() {
		null === this.boundingSphere && (this.boundingSphere = new Sphere);
		const e = this.attributes.position,
			t = this.morphAttributes.position;
		if (e && e.isGLBufferAttribute) return console.error('THREE.BufferGeometry.computeBoundingSphere(): GLBufferAttribute requires a manual bounding sphere. Alternatively set "mesh.frustumCulled" to "false".', this), void this.boundingSphere.set(new Vector3, 1 / 0);
		if (e) {
			const i = this.boundingSphere.center;
			if (_box$2.setFromBufferAttribute(e), t)
				for (let e = 0, i = t.length; e < i; e++) {
					const i = t[e];
					_boxMorphTargets.setFromBufferAttribute(i), this.morphTargetsRelative ? (_vector$4.addVectors(_box$2.min, _boxMorphTargets.min), _box$2.expandByPoint(_vector$4), _vector$4.addVectors(_box$2.max, _boxMorphTargets.max), _box$2.expandByPoint(_vector$4)) : (_box$2.expandByPoint(_boxMorphTargets.min), _box$2.expandByPoint(_boxMorphTargets.max))
				}
			_box$2.getCenter(i);
			let n = 0;
			for (let t = 0, r = e.count; t < r; t++) _vector$4.fromBufferAttribute(e, t), n = Math.max(n, i.distanceToSquared(_vector$4));
			if (t)
				for (let r = 0, a = t.length; r < a; r++) {
					const a = t[r],
						s = this.morphTargetsRelative;
					for (let t = 0, r = a.count; t < r; t++) _vector$4.fromBufferAttribute(a, t), s && (_offset.fromBufferAttribute(e, t), _vector$4.add(_offset)), n = Math.max(n, i.distanceToSquared(_vector$4))
				}
			this.boundingSphere.radius = Math.sqrt(n), isNaN(this.boundingSphere.radius) && console.error('THREE.BufferGeometry.computeBoundingSphere(): Computed radius is NaN. The "position" attribute is likely to have NaN values.', this)
		}
	},
	computeFaceNormals: function() {},
	computeTangents: function() {
		const e = this.index,
			t = this.attributes;
		if (null === e || void 0 === t.position || void 0 === t.normal || void 0 === t.uv) return void console.error("THREE.BufferGeometry: .computeTangents() failed. Missing required attributes (index, position, normal or uv)");
		const i = e.array,
			n = t.position.array,
			r = t.normal.array,
			a = t.uv.array,
			s = n.length / 3;
		void 0 === t.tangent && this.setAttribute("tangent", new BufferAttribute(new Float32Array(4 * s), 4));
		const o = t.tangent.array,
			l = [],
			c = [];
		for (let e = 0; e < s; e++) l[e] = new Vector3, c[e] = new Vector3;
		const h = new Vector3,
			u = new Vector3,
			d = new Vector3,
			p = new Vector2,
			m = new Vector2,
			A = new Vector2,
			g = new Vector3,
			f = new Vector3;

		function v(e, t, i) {
			h.fromArray(n, 3 * e), u.fromArray(n, 3 * t), d.fromArray(n, 3 * i), p.fromArray(a, 2 * e), m.fromArray(a, 2 * t), A.fromArray(a, 2 * i), u.sub(h), d.sub(h), m.sub(p), A.sub(p);
			const r = 1 / (m.x * A.y - A.x * m.y);
			isFinite(r) && (g.copy(u).multiplyScalar(A.y).addScaledVector(d, -m.y).multiplyScalar(r), f.copy(d).multiplyScalar(m.x).addScaledVector(u, -A.x).multiplyScalar(r), l[e].add(g), l[t].add(g), l[i].add(g), c[e].add(f), c[t].add(f), c[i].add(f))
		}
		let y = this.groups;
		0 === y.length && (y = [{
			start: 0,
			count: i.length
		}]);
		for (let e = 0, t = y.length; e < t; ++e) {
			const t = y[e],
				n = t.start;
			for (let e = n, r = n + t.count; e < r; e += 3) v(i[e + 0], i[e + 1], i[e + 2])
		}
		const E = new Vector3,
			_ = new Vector3,
			b = new Vector3,
			x = new Vector3;

		function w(e) {
			b.fromArray(r, 3 * e), x.copy(b);
			const t = l[e];
			E.copy(t), E.sub(b.multiplyScalar(b.dot(t))).normalize(), _.crossVectors(x, t);
			const i = _.dot(c[e]) < 0 ? -1 : 1;
			o[4 * e] = E.x, o[4 * e + 1] = E.y, o[4 * e + 2] = E.z, o[4 * e + 3] = i
		}
		for (let e = 0, t = y.length; e < t; ++e) {
			const t = y[e],
				n = t.start;
			for (let e = n, r = n + t.count; e < r; e += 3) w(i[e + 0]), w(i[e + 1]), w(i[e + 2])
		}
	},
	computeVertexNormals: function() {
		const e = this.index,
			t = this.getAttribute("position");
		if (void 0 !== t) {
			let i = this.getAttribute("normal");
			if (void 0 === i) i = new BufferAttribute(new Float32Array(3 * t.count), 3), this.setAttribute("normal", i);
			else
				for (let e = 0, t = i.count; e < t; e++) i.setXYZ(e, 0, 0, 0);
			const n = new Vector3,
				r = new Vector3,
				a = new Vector3,
				s = new Vector3,
				o = new Vector3,
				l = new Vector3,
				c = new Vector3,
				h = new Vector3;
			if (e)
				for (let u = 0, d = e.count; u < d; u += 3) {
					const d = e.getX(u + 0),
						p = e.getX(u + 1),
						m = e.getX(u + 2);
					n.fromBufferAttribute(t, d), r.fromBufferAttribute(t, p), a.fromBufferAttribute(t, m), c.subVectors(a, r), h.subVectors(n, r), c.cross(h), s.fromBufferAttribute(i, d), o.fromBufferAttribute(i, p), l.fromBufferAttribute(i, m), s.add(c), o.add(c), l.add(c), i.setXYZ(d, s.x, s.y, s.z), i.setXYZ(p, o.x, o.y, o.z), i.setXYZ(m, l.x, l.y, l.z)
				} else
					for (let e = 0, s = t.count; e < s; e += 3) n.fromBufferAttribute(t, e + 0), r.fromBufferAttribute(t, e + 1), a.fromBufferAttribute(t, e + 2), c.subVectors(a, r), h.subVectors(n, r), c.cross(h), i.setXYZ(e + 0, c.x, c.y, c.z), i.setXYZ(e + 1, c.x, c.y, c.z), i.setXYZ(e + 2, c.x, c.y, c.z);
			this.normalizeNormals(), i.needsUpdate = !0
		}
	},
	merge: function(e, t) {
		if (!e || !e.isBufferGeometry) return void console.error("THREE.BufferGeometry.merge(): geometry not an instance of THREE.BufferGeometry.", e);
		void 0 === t && (t = 0, console.warn("THREE.BufferGeometry.merge(): Overwriting original geometry, starting at offset=0. Use BufferGeometryUtils.mergeBufferGeometries() for lossless merge."));
		const i = this.attributes;
		for (const n in i) {
			if (void 0 === e.attributes[n]) continue;
			const r = i[n].array,
				a = e.attributes[n],
				s = a.array,
				o = a.itemSize * t,
				l = Math.min(s.length, r.length - o);
			for (let e = 0, t = o; e < l; e++, t++) r[t] = s[e]
		}
		return this
	},
	normalizeNormals: function() {
		const e = this.attributes.normal;
		for (let t = 0, i = e.count; t < i; t++) _vector$4.fromBufferAttribute(e, t), _vector$4.normalize(), e.setXYZ(t, _vector$4.x, _vector$4.y, _vector$4.z)
	},
	toNonIndexed: function() {
		function e(e, t) {
			const i = e.array,
				n = e.itemSize,
				r = e.normalized,
				a = new i.constructor(t.length * n);
			let s = 0,
				o = 0;
			for (let e = 0, r = t.length; e < r; e++) {
				s = t[e] * n;
				for (let e = 0; e < n; e++) a[o++] = i[s++]
			}
			return new BufferAttribute(a, n, r)
		}
		if (null === this.index) return console.warn("THREE.BufferGeometry.toNonIndexed(): BufferGeometry is already non-indexed."), this;
		const t = new BufferGeometry,
			i = this.index.array,
			n = this.attributes;
		for (const r in n) {
			const a = e(n[r], i);
			t.setAttribute(r, a)
		}
		const r = this.morphAttributes;
		for (const n in r) {
			const a = [],
				s = r[n];
			for (let t = 0, n = s.length; t < n; t++) {
				const n = e(s[t], i);
				a.push(n)
			}
			t.morphAttributes[n] = a
		}
		t.morphTargetsRelative = this.morphTargetsRelative;
		const a = this.groups;
		for (let e = 0, i = a.length; e < i; e++) {
			const i = a[e];
			t.addGroup(i.start, i.count, i.materialIndex)
		}
		return t
	},
	toJSON: function() {
		const e = {
			metadata: {
				version: 4.5,
				type: "BufferGeometry",
				generator: "BufferGeometry.toJSON"
			}
		};
		if (e.uuid = this.uuid, e.type = this.type, "" !== this.name && (e.name = this.name), Object.keys(this.userData).length > 0 && (e.userData = this.userData), void 0 !== this.parameters) {
			const t = this.parameters;
			for (const i in t) void 0 !== t[i] && (e[i] = t[i]);
			return e
		}
		e.data = {
			attributes: {}
		};
		const t = this.index;
		null !== t && (e.data.index = {
			type: t.array.constructor.name,
			array: Array.prototype.slice.call(t.array)
		});
		const i = this.attributes;
		for (const t in i) {
			const n = i[t],
				r = n.toJSON(e.data);
			"" !== n.name && (r.name = n.name), e.data.attributes[t] = r
		}
		const n = {};
		let r = !1;
		for (const t in this.morphAttributes) {
			const i = this.morphAttributes[t],
				a = [];
			for (let t = 0, n = i.length; t < n; t++) {
				const n = i[t],
					r = n.toJSON(e.data);
				"" !== n.name && (r.name = n.name), a.push(r)
			}
			a.length > 0 && (n[t] = a, r = !0)
		}
		r && (e.data.morphAttributes = n, e.data.morphTargetsRelative = this.morphTargetsRelative);
		const a = this.groups;
		a.length > 0 && (e.data.groups = JSON.parse(JSON.stringify(a)));
		const s = this.boundingSphere;
		return null !== s && (e.data.boundingSphere = {
			center: s.center.toArray(),
			radius: s.radius
		}), e
	},
	clone: function() {
		return (new BufferGeometry).copy(this)
	},
	copy: function(e) {
		this.index = null, this.attributes = {}, this.morphAttributes = {}, this.groups = [], this.boundingBox = null, this.boundingSphere = null;
		const t = {};
		this.name = e.name;
		const i = e.index;
		null !== i && this.setIndex(i.clone(t));
		const n = e.attributes;
		for (const e in n) {
			const i = n[e];
			this.setAttribute(e, i.clone(t))
		}
		const r = e.morphAttributes;
		for (const e in r) {
			const i = [],
				n = r[e];
			for (let e = 0, r = n.length; e < r; e++) i.push(n[e].clone(t));
			this.morphAttributes[e] = i
		}
		this.morphTargetsRelative = e.morphTargetsRelative;
		const a = e.groups;
		for (let e = 0, t = a.length; e < t; e++) {
			const t = a[e];
			this.addGroup(t.start, t.count, t.materialIndex)
		}
		const s = e.boundingBox;
		null !== s && (this.boundingBox = s.clone());
		const o = e.boundingSphere;
		return null !== o && (this.boundingSphere = o.clone()), this.drawRange.start = e.drawRange.start, this.drawRange.count = e.drawRange.count, this.userData = e.userData, this
	},
	dispose: function() {
		this.dispatchEvent({
			type: "dispose"
		})
	}
});
const _inverseMatrix = new Matrix4,
	_ray = new Ray,
	_sphere = new Sphere,
	_vA = new Vector3,
	_vB = new Vector3,
	_vC = new Vector3,
	_tempA = new Vector3,
	_tempB = new Vector3,
	_tempC = new Vector3,
	_morphA = new Vector3,
	_morphB = new Vector3,
	_morphC = new Vector3,
	_uvA = new Vector2,
	_uvB = new Vector2,
	_uvC = new Vector2,
	_intersectionPoint = new Vector3,
	_intersectionPointWorld = new Vector3;

function Mesh(e = new BufferGeometry, t = new MeshBasicMaterial) {
	Object3D.call(this), this.type = "Mesh", this.geometry = e, this.material = t, this.updateMorphTargets()
}

function checkIntersection(e, t, i, n, r, a, s, o) {
	let l;
	if (l = 1 === t.side ? n.intersectTriangle(s, a, r, !0, o) : n.intersectTriangle(r, a, s, 2 !== t.side, o), null === l) return null;
	_intersectionPointWorld.copy(o), _intersectionPointWorld.applyMatrix4(e.matrixWorld);
	const c = i.ray.origin.distanceTo(_intersectionPointWorld);
	return c < i.near || c > i.far ? null : {
		distance: c,
		point: _intersectionPointWorld.clone(),
		object: e
	}
}

function checkBufferGeometryIntersection(e, t, i, n, r, a, s, o, l, c, h, u) {
	_vA.fromBufferAttribute(r, c), _vB.fromBufferAttribute(r, h), _vC.fromBufferAttribute(r, u);
	const d = e.morphTargetInfluences;
	if (t.morphTargets && a && d) {
		_morphA.set(0, 0, 0), _morphB.set(0, 0, 0), _morphC.set(0, 0, 0);
		for (let e = 0, t = a.length; e < t; e++) {
			const t = d[e],
				i = a[e];
			0 !== t && (_tempA.fromBufferAttribute(i, c), _tempB.fromBufferAttribute(i, h), _tempC.fromBufferAttribute(i, u), s ? (_morphA.addScaledVector(_tempA, t), _morphB.addScaledVector(_tempB, t), _morphC.addScaledVector(_tempC, t)) : (_morphA.addScaledVector(_tempA.sub(_vA), t), _morphB.addScaledVector(_tempB.sub(_vB), t), _morphC.addScaledVector(_tempC.sub(_vC), t)))
		}
		_vA.add(_morphA), _vB.add(_morphB), _vC.add(_morphC)
	}
	e.isSkinnedMesh && t.skinning && (e.boneTransform(c, _vA), e.boneTransform(h, _vB), e.boneTransform(u, _vC));
	const p = checkIntersection(e, t, i, n, _vA, _vB, _vC, _intersectionPoint);
	if (p) {
		o && (_uvA.fromBufferAttribute(o, c), _uvB.fromBufferAttribute(o, h), _uvC.fromBufferAttribute(o, u), p.uv = Triangle.getUV(_intersectionPoint, _vA, _vB, _vC, _uvA, _uvB, _uvC, new Vector2)), l && (_uvA.fromBufferAttribute(l, c), _uvB.fromBufferAttribute(l, h), _uvC.fromBufferAttribute(l, u), p.uv2 = Triangle.getUV(_intersectionPoint, _vA, _vB, _vC, _uvA, _uvB, _uvC, new Vector2));
		const e = {
			a: c,
			b: c,
			c: u,
			normal: new Vector3,
			materialIndex: 0
		};
		Triangle.getNormal(_vA, _vB, _vC, e.normal), p.face = e
	}
	return p
}
Mesh.prototype = Object.assign(Object.create(Object3D.prototype), {
	constructor: Mesh,
	isMesh: !0,
	copy: function(e) {
		return Object3D.prototype.copy.call(this, e), void 0 !== e.morphTargetInfluences && (this.morphTargetInfluences = e.morphTargetInfluences.slice()), void 0 !== e.morphTargetDictionary && (this.morphTargetDictionary = Object.assign({}, e.morphTargetDictionary)), this.material = e.material, this.geometry = e.geometry, this
	},
	updateMorphTargets: function() {
		const e = this.geometry;
		if (e.isBufferGeometry) {
			const t = e.morphAttributes,
				i = Object.keys(t);
			if (i.length > 0) {
				const e = t[i[0]];
				if (void 0 !== e) {
					this.morphTargetInfluences = [], this.morphTargetDictionary = {};
					for (let t = 0, i = e.length; t < i; t++) {
						const i = e[t].name || String(t);
						this.morphTargetInfluences.push(0), this.morphTargetDictionary[i] = t
					}
				}
			}
		} else {
			const t = e.morphTargets;
			void 0 !== t && t.length > 0 && console.error("THREE.Mesh.updateMorphTargets() no longer supports THREE.Geometry. Use THREE.BufferGeometry instead.")
		}
	},
	raycast: function(e, t) {
		const i = this.geometry,
			n = this.material,
			r = this.matrixWorld;
		if (void 0 === n) return;
		if (null === i.boundingSphere && i.computeBoundingSphere(), _sphere.copy(i.boundingSphere), _sphere.applyMatrix4(r), !1 === e.ray.intersectsSphere(_sphere)) return;
		if (_inverseMatrix.copy(r).invert(), _ray.copy(e.ray).applyMatrix4(_inverseMatrix), null !== i.boundingBox && !1 === _ray.intersectsBox(i.boundingBox)) return;
		let a;
		if (i.isBufferGeometry) {
			const r = i.index,
				s = i.attributes.position,
				o = i.morphAttributes.position,
				l = i.morphTargetsRelative,
				c = i.attributes.uv,
				h = i.attributes.uv2,
				u = i.groups,
				d = i.drawRange;
			if (null !== r)
				if (Array.isArray(n))
					for (let i = 0, p = u.length; i < p; i++) {
						const p = u[i],
							m = n[p.materialIndex];
						for (let i = Math.max(p.start, d.start), n = Math.min(p.start + p.count, d.start + d.count); i < n; i += 3) {
							const n = r.getX(i),
								u = r.getX(i + 1),
								d = r.getX(i + 2);
							a = checkBufferGeometryIntersection(this, m, e, _ray, s, o, l, c, h, n, u, d), a && (a.faceIndex = Math.floor(i / 3), a.face.materialIndex = p.materialIndex, t.push(a))
						}
					} else {
						for (let i = Math.max(0, d.start), u = Math.min(r.count, d.start + d.count); i < u; i += 3) {
							const u = r.getX(i),
								d = r.getX(i + 1),
								p = r.getX(i + 2);
							a = checkBufferGeometryIntersection(this, n, e, _ray, s, o, l, c, h, u, d, p), a && (a.faceIndex = Math.floor(i / 3), t.push(a))
						}
					} else if (void 0 !== s)
						if (Array.isArray(n))
							for (let i = 0, r = u.length; i < r; i++) {
								const r = u[i],
									p = n[r.materialIndex];
								for (let i = Math.max(r.start, d.start), n = Math.min(r.start + r.count, d.start + d.count); i < n; i += 3) {
									a = checkBufferGeometryIntersection(this, p, e, _ray, s, o, l, c, h, i, i + 1, i + 2), a && (a.faceIndex = Math.floor(i / 3), a.face.materialIndex = r.materialIndex, t.push(a))
								}
							} else {
								for (let i = Math.max(0, d.start), r = Math.min(s.count, d.start + d.count); i < r; i += 3) {
									a = checkBufferGeometryIntersection(this, n, e, _ray, s, o, l, c, h, i, i + 1, i + 2), a && (a.faceIndex = Math.floor(i / 3), t.push(a))
								}
							}
		} else i.isGeometry && console.error("THREE.Mesh.raycast() no longer supports THREE.Geometry. Use THREE.BufferGeometry instead.")
	}
});
class BoxGeometry extends BufferGeometry {
	constructor(e = 1, t = 1, i = 1, n = 1, r = 1, a = 1) {
		super(), this.type = "BoxGeometry", this.parameters = {
			width: e,
			height: t,
			depth: i,
			widthSegments: n,
			heightSegments: r,
			depthSegments: a
		};
		const s = this;
		n = Math.floor(n), r = Math.floor(r), a = Math.floor(a);
		const o = [],
			l = [],
			c = [],
			h = [];
		let u = 0,
			d = 0;

		function p(e, t, i, n, r, a, p, m, A, g, f) {
			const v = a / A,
				y = p / g,
				E = a / 2,
				_ = p / 2,
				b = m / 2,
				x = A + 1,
				w = g + 1;
			let C = 0,
				S = 0;
			const I = new Vector3;
			for (let a = 0; a < w; a++) {
				const s = a * y - _;
				for (let o = 0; o < x; o++) {
					const u = o * v - E;
					I[e] = u * n, I[t] = s * r, I[i] = b, l.push(I.x, I.y, I.z), I[e] = 0, I[t] = 0, I[i] = m > 0 ? 1 : -1, c.push(I.x, I.y, I.z), h.push(o / A), h.push(1 - a / g), C += 1
				}
			}
			for (let e = 0; e < g; e++)
				for (let t = 0; t < A; t++) {
					const i = u + t + x * e,
						n = u + t + x * (e + 1),
						r = u + (t + 1) + x * (e + 1),
						a = u + (t + 1) + x * e;
					o.push(i, n, a), o.push(n, r, a), S += 6
				}
			s.addGroup(d, S, f), d += S, u += C
		}
		p("z", "y", "x", -1, -1, i, t, e, a, r, 0), p("z", "y", "x", 1, -1, i, t, -e, a, r, 1), p("x", "z", "y", 1, 1, e, i, t, n, a, 2), p("x", "z", "y", 1, -1, e, i, -t, n, a, 3), p("x", "y", "z", 1, -1, e, t, i, n, r, 4), p("x", "y", "z", -1, -1, e, t, -i, n, r, 5), this.setIndex(o), this.setAttribute("position", new Float32BufferAttribute(l, 3)), this.setAttribute("normal", new Float32BufferAttribute(c, 3)), this.setAttribute("uv", new Float32BufferAttribute(h, 2))
	}
}

function cloneUniforms(e) {
	const t = {};
	for (const i in e) {
		t[i] = {};
		for (const n in e[i]) {
			const r = e[i][n];
			r && (r.isColor || r.isMatrix3 || r.isMatrix4 || r.isVector2 || r.isVector3 || r.isVector4 || r.isTexture || r.isQuaternion) ? t[i][n] = r.clone() : Array.isArray(r) ? t[i][n] = r.slice() : t[i][n] = r
		}
	}
	return t
}

function mergeUniforms(e) {
	const t = {};
	for (let i = 0; i < e.length; i++) {
		const n = cloneUniforms(e[i]);
		for (const e in n) t[e] = n[e]
	}
	return t
}
const UniformsUtils = {
	clone: cloneUniforms,
	merge: mergeUniforms
};
var default_vertex = "void main() {\n\tgl_Position = projectionMatrix * modelViewMatrix * vec4( position, 1.0 );\n}",
	default_fragment = "void main() {\n\tgl_FragColor = vec4( 1.0, 0.0, 0.0, 1.0 );\n}";

function ShaderMaterial(e) {
	Material$1.call(this), this.type = "ShaderMaterial", this.defines = {}, this.uniforms = {}, this.vertexShader = default_vertex, this.fragmentShader = default_fragment, this.linewidth = 1, this.wireframe = !1, this.wireframeLinewidth = 1, this.fog = !1, this.lights = !1, this.clipping = !1, this.skinning = !1, this.morphTargets = !1, this.morphNormals = !1, this.extensions = {
		derivatives: !1,
		fragDepth: !1,
		drawBuffers: !1,
		shaderTextureLOD: !1
	}, this.defaultAttributeValues = {
		color: [1, 1, 1],
		uv: [0, 0],
		uv2: [0, 0]
	}, this.index0AttributeName = void 0, this.uniformsNeedUpdate = !1, this.glslVersion = null, void 0 !== e && (void 0 !== e.attributes && console.error("THREE.ShaderMaterial: attributes should now be defined in THREE.BufferGeometry instead."), this.setValues(e))
}

function Camera() {
	Object3D.call(this), this.type = "Camera", this.matrixWorldInverse = new Matrix4, this.projectionMatrix = new Matrix4, this.projectionMatrixInverse = new Matrix4
}

function PerspectiveCamera(e = 50, t = 1, i = .1, n = 2e3) {
	Camera.call(this), this.type = "PerspectiveCamera", this.fov = e, this.zoom = 1, this.near = i, this.far = n, this.focus = 10, this.aspect = t, this.view = null, this.filmGauge = 35, this.filmOffset = 0, this.updateProjectionMatrix()
}
ShaderMaterial.prototype = Object.create(Material$1.prototype), ShaderMaterial.prototype.constructor = ShaderMaterial, ShaderMaterial.prototype.isShaderMaterial = !0, ShaderMaterial.prototype.copy = function(e) {
	return Material$1.prototype.copy.call(this, e), this.fragmentShader = e.fragmentShader, this.vertexShader = e.vertexShader, this.uniforms = cloneUniforms(e.uniforms), this.defines = Object.assign({}, e.defines), this.wireframe = e.wireframe, this.wireframeLinewidth = e.wireframeLinewidth, this.lights = e.lights, this.clipping = e.clipping, this.skinning = e.skinning, this.morphTargets = e.morphTargets, this.morphNormals = e.morphNormals, this.extensions = Object.assign({}, e.extensions), this.glslVersion = e.glslVersion, this
}, ShaderMaterial.prototype.toJSON = function(e) {
	const t = Material$1.prototype.toJSON.call(this, e);
	t.glslVersion = this.glslVersion, t.uniforms = {};
	for (const i in this.uniforms) {
		const n = this.uniforms[i].value;
		n && n.isTexture ? t.uniforms[i] = {
			type: "t",
			value: n.toJSON(e).uuid
		} : n && n.isColor ? t.uniforms[i] = {
			type: "c",
			value: n.getHex()
		} : n && n.isVector2 ? t.uniforms[i] = {
			type: "v2",
			value: n.toArray()
		} : n && n.isVector3 ? t.uniforms[i] = {
			type: "v3",
			value: n.toArray()
		} : n && n.isVector4 ? t.uniforms[i] = {
			type: "v4",
			value: n.toArray()
		} : n && n.isMatrix3 ? t.uniforms[i] = {
			type: "m3",
			value: n.toArray()
		} : n && n.isMatrix4 ? t.uniforms[i] = {
			type: "m4",
			value: n.toArray()
		} : t.uniforms[i] = {
			value: n
		}
	}
	Object.keys(this.defines).length > 0 && (t.defines = this.defines), t.vertexShader = this.vertexShader, t.fragmentShader = this.fragmentShader;
	const i = {};
	for (const e in this.extensions) !0 === this.extensions[e] && (i[e] = !0);
	return Object.keys(i).length > 0 && (t.extensions = i), t
}, Camera.prototype = Object.assign(Object.create(Object3D.prototype), {
	constructor: Camera,
	isCamera: !0,
	copy: function(e, t) {
		return Object3D.prototype.copy.call(this, e, t), this.matrixWorldInverse.copy(e.matrixWorldInverse), this.projectionMatrix.copy(e.projectionMatrix), this.projectionMatrixInverse.copy(e.projectionMatrixInverse), this
	},
	getWorldDirection: function(e) {
		void 0 === e && (console.warn("THREE.Camera: .getWorldDirection() target is now required"), e = new Vector3), this.updateWorldMatrix(!0, !1);
		const t = this.matrixWorld.elements;
		return e.set(-t[8], -t[9], -t[10]).normalize()
	},
	updateMatrixWorld: function(e) {
		Object3D.prototype.updateMatrixWorld.call(this, e), this.matrixWorldInverse.copy(this.matrixWorld).invert()
	},
	updateWorldMatrix: function(e, t) {
		Object3D.prototype.updateWorldMatrix.call(this, e, t), this.matrixWorldInverse.copy(this.matrixWorld).invert()
	},
	clone: function() {
		return (new this.constructor).copy(this)
	}
}), PerspectiveCamera.prototype = Object.assign(Object.create(Camera.prototype), {
	constructor: PerspectiveCamera,
	isPerspectiveCamera: !0,
	copy: function(e, t) {
		return Camera.prototype.copy.call(this, e, t), this.fov = e.fov, this.zoom = e.zoom, this.near = e.near, this.far = e.far, this.focus = e.focus, this.aspect = e.aspect, this.view = null === e.view ? null : Object.assign({}, e.view), this.filmGauge = e.filmGauge, this.filmOffset = e.filmOffset, this
	},
	setFocalLength: function(e) {
		const t = .5 * this.getFilmHeight() / e;
		this.fov = 2 * MathUtils.RAD2DEG * Math.atan(t), this.updateProjectionMatrix()
	},
	getFocalLength: function() {
		const e = Math.tan(.5 * MathUtils.DEG2RAD * this.fov);
		return .5 * this.getFilmHeight() / e
	},
	getEffectiveFOV: function() {
		return 2 * MathUtils.RAD2DEG * Math.atan(Math.tan(.5 * MathUtils.DEG2RAD * this.fov) / this.zoom)
	},
	getFilmWidth: function() {
		return this.filmGauge * Math.min(this.aspect, 1)
	},
	getFilmHeight: function() {
		return this.filmGauge / Math.max(this.aspect, 1)
	},
	setViewOffset: function(e, t, i, n, r, a) {
		this.aspect = e / t, null === this.view && (this.view = {
			enabled: !0,
			fullWidth: 1,
			fullHeight: 1,
			offsetX: 0,
			offsetY: 0,
			width: 1,
			height: 1
		}), this.view.enabled = !0, this.view.fullWidth = e, this.view.fullHeight = t, this.view.offsetX = i, this.view.offsetY = n, this.view.width = r, this.view.height = a, this.updateProjectionMatrix()
	},
	clearViewOffset: function() {
		null !== this.view && (this.view.enabled = !1), this.updateProjectionMatrix()
	},
	updateProjectionMatrix: function() {
		const e = this.near;
		let t = e * Math.tan(.5 * MathUtils.DEG2RAD * this.fov) / this.zoom,
			i = 2 * t,
			n = this.aspect * i,
			r = -.5 * n;
		const a = this.view;
		if (null !== this.view && this.view.enabled) {
			const e = a.fullWidth,
				s = a.fullHeight;
			r += a.offsetX * n / e, t -= a.offsetY * i / s, n *= a.width / e, i *= a.height / s
		}
		const s = this.filmOffset;
		0 !== s && (r += e * s / this.getFilmWidth()), this.projectionMatrix.makePerspective(r, r + n, t, t - i, e, this.far), this.projectionMatrixInverse.copy(this.projectionMatrix).invert()
	},
	toJSON: function(e) {
		const t = Object3D.prototype.toJSON.call(this, e);
		return t.object.fov = this.fov, t.object.zoom = this.zoom, t.object.near = this.near, t.object.far = this.far, t.object.focus = this.focus, t.object.aspect = this.aspect, null !== this.view && (t.object.view = Object.assign({}, this.view)), t.object.filmGauge = this.filmGauge, t.object.filmOffset = this.filmOffset, t
	}
});
const fov = 90,
	aspect = 1;
class CubeCamera extends Object3D {
	constructor(e, t, i) {
		if (super(), this.type = "CubeCamera", !0 !== i.isWebGLCubeRenderTarget) return void console.error("THREE.CubeCamera: The constructor now expects an instance of WebGLCubeRenderTarget as third parameter.");
		this.renderTarget = i;
		const n = new PerspectiveCamera(90, 1, e, t);
		n.layers = this.layers, n.up.set(0, -1, 0), n.lookAt(new Vector3(1, 0, 0)), this.add(n);
		const r = new PerspectiveCamera(90, 1, e, t);
		r.layers = this.layers, r.up.set(0, -1, 0), r.lookAt(new Vector3(-1, 0, 0)), this.add(r);
		const a = new PerspectiveCamera(90, 1, e, t);
		a.layers = this.layers, a.up.set(0, 0, 1), a.lookAt(new Vector3(0, 1, 0)), this.add(a);
		const s = new PerspectiveCamera(90, 1, e, t);
		s.layers = this.layers, s.up.set(0, 0, -1), s.lookAt(new Vector3(0, -1, 0)), this.add(s);
		const o = new PerspectiveCamera(90, 1, e, t);
		o.layers = this.layers, o.up.set(0, -1, 0), o.lookAt(new Vector3(0, 0, 1)), this.add(o);
		const l = new PerspectiveCamera(90, 1, e, t);
		l.layers = this.layers, l.up.set(0, -1, 0), l.lookAt(new Vector3(0, 0, -1)), this.add(l)
	}
	update(e, t) {
		null === this.parent && this.updateMatrixWorld();
		const i = this.renderTarget,
			[n, r, a, s, o, l] = this.children,
			c = e.xr.enabled,
			h = e.getRenderTarget();
		e.xr.enabled = !1;
		const u = i.texture.generateMipmaps;
		i.texture.generateMipmaps = !1, e.setRenderTarget(i, 0), e.render(t, n), e.setRenderTarget(i, 1), e.render(t, r), e.setRenderTarget(i, 2), e.render(t, a), e.setRenderTarget(i, 3), e.render(t, s), e.setRenderTarget(i, 4), e.render(t, o), i.texture.generateMipmaps = u, e.setRenderTarget(i, 5), e.render(t, l), e.setRenderTarget(h), e.xr.enabled = c
	}
}
class CubeTexture extends Texture$1 {
	constructor(e, t, i, n, r, a, s, o, l, c) {
		super(e = void 0 !== e ? e : [], t = void 0 !== t ? t : 301, i, n, r, a, s = void 0 !== s ? s : 1022, o, l, c), this._needsFlipEnvMap = !0, this.flipY = !1
	}
	get images() {
		return this.image
	}
	set images(e) {
		this.image = e
	}
}
CubeTexture.prototype.isCubeTexture = !0;
class WebGLCubeRenderTarget extends WebGLRenderTarget {
	constructor(e, t, i) {
		Number.isInteger(t) && (console.warn("THREE.WebGLCubeRenderTarget: constructor signature is now WebGLCubeRenderTarget( size, options )"), t = i), super(e, e, t), t = t || {}, this.texture = new CubeTexture(void 0, t.mapping, t.wrapS, t.wrapT, t.magFilter, t.minFilter, t.format, t.type, t.anisotropy, t.encoding), this.texture.generateMipmaps = void 0 !== t.generateMipmaps && t.generateMipmaps, this.texture.minFilter = void 0 !== t.minFilter ? t.minFilter : 1006, this.texture._needsFlipEnvMap = !1
	}
	fromEquirectangularTexture(e, t) {
		this.texture.type = t.type, this.texture.format = 1023, this.texture.encoding = t.encoding, this.texture.generateMipmaps = t.generateMipmaps, this.texture.minFilter = t.minFilter, this.texture.magFilter = t.magFilter;
		const i = {
				uniforms: {
					tEquirect: {
						value: null
					}
				},
				vertexShader: "\n\n\t\t\t\tvarying vec3 vWorldDirection;\n\n\t\t\t\tvec3 transformDirection( in vec3 dir, in mat4 matrix ) {\n\n\t\t\t\t\treturn normalize( ( matrix * vec4( dir, 0.0 ) ).xyz );\n\n\t\t\t\t}\n\n\t\t\t\tvoid main() {\n\n\t\t\t\t\tvWorldDirection = transformDirection( position, modelMatrix );\n\n\t\t\t\t\t#include <begin_vertex>\n\t\t\t\t\t#include <project_vertex>\n\n\t\t\t\t}\n\t\t\t",
				fragmentShader: "\n\n\t\t\t\tuniform sampler2D tEquirect;\n\n\t\t\t\tvarying vec3 vWorldDirection;\n\n\t\t\t\t#include <common>\n\n\t\t\t\tvoid main() {\n\n\t\t\t\t\tvec3 direction = normalize( vWorldDirection );\n\n\t\t\t\t\tvec2 sampleUV = equirectUv( direction );\n\n\t\t\t\t\tgl_FragColor = texture2D( tEquirect, sampleUV );\n\n\t\t\t\t}\n\t\t\t"
			},
			n = new BoxGeometry(5, 5, 5),
			r = new ShaderMaterial({
				name: "CubemapFromEquirect",
				uniforms: cloneUniforms(i.uniforms),
				vertexShader: i.vertexShader,
				fragmentShader: i.fragmentShader,
				side: 1,
				blending: 0
			});
		r.uniforms.tEquirect.value = t;
		const a = new Mesh(n, r),
			s = t.minFilter;
		1008 === t.minFilter && (t.minFilter = 1006);
		return new CubeCamera(1, 10, this).update(e, a), t.minFilter = s, a.geometry.dispose(), a.material.dispose(), this
	}
	clear(e, t, i, n) {
		const r = e.getRenderTarget();
		for (let r = 0; r < 6; r++) e.setRenderTarget(this, r), e.clear(t, i, n);
		e.setRenderTarget(r)
	}
}
WebGLCubeRenderTarget.prototype.isWebGLCubeRenderTarget = !0;
class DataTexture extends Texture$1 {
	constructor(e, t, i, n, r, a, s, o, l, c, h, u) {
		super(null, a, s, o, l, c, n, r, h, u), this.image = {
			data: e || null,
			width: t || 1,
			height: i || 1
		}, this.magFilter = void 0 !== l ? l : 1003, this.minFilter = void 0 !== c ? c : 1003, this.generateMipmaps = !1, this.flipY = !1, this.unpackAlignment = 1, this.needsUpdate = !0
	}
}
DataTexture.prototype.isDataTexture = !0;
const _sphere$1 = new Sphere,
	_vector$5 = new Vector3;
class Frustum {
	constructor(e = new Plane, t = new Plane, i = new Plane, n = new Plane, r = new Plane, a = new Plane) {
		this.planes = [e, t, i, n, r, a]
	}
	set(e, t, i, n, r, a) {
		const s = this.planes;
		return s[0].copy(e), s[1].copy(t), s[2].copy(i), s[3].copy(n), s[4].copy(r), s[5].copy(a), this
	}
	copy(e) {
		const t = this.planes;
		for (let i = 0; i < 6; i++) t[i].copy(e.planes[i]);
		return this
	}
	setFromProjectionMatrix(e) {
		const t = this.planes,
			i = e.elements,
			n = i[0],
			r = i[1],
			a = i[2],
			s = i[3],
			o = i[4],
			l = i[5],
			c = i[6],
			h = i[7],
			u = i[8],
			d = i[9],
			p = i[10],
			m = i[11],
			A = i[12],
			g = i[13],
			f = i[14],
			v = i[15];
		return t[0].setComponents(s - n, h - o, m - u, v - A).normalize(), t[1].setComponents(s + n, h + o, m + u, v + A).normalize(), t[2].setComponents(s + r, h + l, m + d, v + g).normalize(), t[3].setComponents(s - r, h - l, m - d, v - g).normalize(), t[4].setComponents(s - a, h - c, m - p, v - f).normalize(), t[5].setComponents(s + a, h + c, m + p, v + f).normalize(), this
	}
	intersectsObject(e) {
		const t = e.geometry;
		return null === t.boundingSphere && t.computeBoundingSphere(), _sphere$1.copy(t.boundingSphere).applyMatrix4(e.matrixWorld), this.intersectsSphere(_sphere$1)
	}
	intersectsSprite(e) {
		return _sphere$1.center.set(0, 0, 0), _sphere$1.radius = .7071067811865476, _sphere$1.applyMatrix4(e.matrixWorld), this.intersectsSphere(_sphere$1)
	}
	intersectsSphere(e) {
		const t = this.planes,
			i = e.center,
			n = -e.radius;
		for (let e = 0; e < 6; e++) {
			if (t[e].distanceToPoint(i) < n) return !1
		}
		return !0
	}
	intersectsBox(e) {
		const t = this.planes;
		for (let i = 0; i < 6; i++) {
			const n = t[i];
			if (_vector$5.x = n.normal.x > 0 ? e.max.x : e.min.x, _vector$5.y = n.normal.y > 0 ? e.max.y : e.min.y, _vector$5.z = n.normal.z > 0 ? e.max.z : e.min.z, n.distanceToPoint(_vector$5) < 0) return !1
		}
		return !0
	}
	containsPoint(e) {
		const t = this.planes;
		for (let i = 0; i < 6; i++)
			if (t[i].distanceToPoint(e) < 0) return !1;
		return !0
	}
	clone() {
		return (new this.constructor).copy(this)
	}
}

function WebGLAnimation() {
	let e = null,
		t = !1,
		i = null,
		n = null;

	function r(t, a) {
		i(t, a), n = e.requestAnimationFrame(r)
	}
	return {
		start: function() {
			!0 !== t && null !== i && (n = e.requestAnimationFrame(r), t = !0)
		},
		stop: function() {
			e.cancelAnimationFrame(n), t = !1
		},
		setAnimationLoop: function(e) {
			i = e
		},
		setContext: function(t) {
			e = t
		}
	}
}

function WebGLAttributes(e, t) {
	const i = t.isWebGL2,
		n = new WeakMap;
	return {
		get: function(e) {
			return e.isInterleavedBufferAttribute && (e = e.data), n.get(e)
		},
		remove: function(t) {
			t.isInterleavedBufferAttribute && (t = t.data);
			const i = n.get(t);
			i && (e.deleteBuffer(i.buffer), n.delete(t))
		},
		update: function(t, r) {
			if (t.isGLBufferAttribute) {
				const e = n.get(t);
				return void((!e || e.version < t.version) && n.set(t, {
					buffer: t.buffer,
					type: t.type,
					bytesPerElement: t.elementSize,
					version: t.version
				}))
			}
			t.isInterleavedBufferAttribute && (t = t.data);
			const a = n.get(t);
			void 0 === a ? n.set(t, function(t, n) {
				const r = t.array,
					a = t.usage,
					s = e.createBuffer();
				e.bindBuffer(n, s), e.bufferData(n, r, a), t.onUploadCallback();
				let o = 5126;
				return r instanceof Float32Array ? o = 5126 : r instanceof Float64Array ? console.warn("THREE.WebGLAttributes: Unsupported data buffer format: Float64Array.") : r instanceof Uint16Array ? t.isFloat16BufferAttribute ? i ? o = 5131 : console.warn("THREE.WebGLAttributes: Usage of Float16BufferAttribute requires WebGL2.") : o = 5123 : r instanceof Int16Array ? o = 5122 : r instanceof Uint32Array ? o = 5125 : r instanceof Int32Array ? o = 5124 : r instanceof Int8Array ? o = 5120 : r instanceof Uint8Array && (o = 5121), {
					buffer: s,
					type: o,
					bytesPerElement: r.BYTES_PER_ELEMENT,
					version: t.version
				}
			}(t, r)) : a.version < t.version && (! function(t, n, r) {
				const a = n.array,
					s = n.updateRange;
				e.bindBuffer(r, t), -1 === s.count ? e.bufferSubData(r, 0, a) : (i ? e.bufferSubData(r, s.offset * a.BYTES_PER_ELEMENT, a, s.offset, s.count) : e.bufferSubData(r, s.offset * a.BYTES_PER_ELEMENT, a.subarray(s.offset, s.offset + s.count)), s.count = -1)
			}(a.buffer, t, r), a.version = t.version)
		}
	}
}
class PlaneGeometry extends BufferGeometry {
	constructor(e = 1, t = 1, i = 1, n = 1) {
		super(), this.type = "PlaneGeometry", this.parameters = {
			width: e,
			height: t,
			widthSegments: i,
			heightSegments: n
		};
		const r = e / 2,
			a = t / 2,
			s = Math.floor(i),
			o = Math.floor(n),
			l = s + 1,
			c = o + 1,
			h = e / s,
			u = t / o,
			d = [],
			p = [],
			m = [],
			A = [];
		for (let e = 0; e < c; e++) {
			const t = e * u - a;
			for (let i = 0; i < l; i++) {
				const n = i * h - r;
				p.push(n, -t, 0), m.push(0, 0, 1), A.push(i / s), A.push(1 - e / o)
			}
		}
		for (let e = 0; e < o; e++)
			for (let t = 0; t < s; t++) {
				const i = t + l * e,
					n = t + l * (e + 1),
					r = t + 1 + l * (e + 1),
					a = t + 1 + l * e;
				d.push(i, n, a), d.push(n, r, a)
			}
		this.setIndex(d), this.setAttribute("position", new Float32BufferAttribute(p, 3)), this.setAttribute("normal", new Float32BufferAttribute(m, 3)), this.setAttribute("uv", new Float32BufferAttribute(A, 2))
	}
}
var alphamap_fragment = "#ifdef USE_ALPHAMAP\n\tdiffuseColor.a *= texture2D( alphaMap, vUv ).g;\n#endif",
	alphamap_pars_fragment = "#ifdef USE_ALPHAMAP\n\tuniform sampler2D alphaMap;\n#endif",
	alphatest_fragment = "#ifdef ALPHATEST\n\tif ( diffuseColor.a < ALPHATEST ) discard;\n#endif",
	aomap_fragment = "#ifdef USE_AOMAP\n\tfloat ambientOcclusion = ( texture2D( aoMap, vUv2 ).r - 1.0 ) * aoMapIntensity + 1.0;\n\treflectedLight.indirectDiffuse *= ambientOcclusion;\n\t#if defined( USE_ENVMAP ) && defined( STANDARD )\n\t\tfloat dotNV = saturate( dot( geometry.normal, geometry.viewDir ) );\n\t\treflectedLight.indirectSpecular *= computeSpecularOcclusion( dotNV, ambientOcclusion, material.specularRoughness );\n\t#endif\n#endif",
	aomap_pars_fragment = "#ifdef USE_AOMAP\n\tuniform sampler2D aoMap;\n\tuniform float aoMapIntensity;\n#endif",
	begin_vertex = "vec3 transformed = vec3( position );",
	beginnormal_vertex = "vec3 objectNormal = vec3( normal );\n#ifdef USE_TANGENT\n\tvec3 objectTangent = vec3( tangent.xyz );\n#endif",
	bsdfs = "vec2 integrateSpecularBRDF( const in float dotNV, const in float roughness ) {\n\tconst vec4 c0 = vec4( - 1, - 0.0275, - 0.572, 0.022 );\n\tconst vec4 c1 = vec4( 1, 0.0425, 1.04, - 0.04 );\n\tvec4 r = roughness * c0 + c1;\n\tfloat a004 = min( r.x * r.x, exp2( - 9.28 * dotNV ) ) * r.x + r.y;\n\treturn vec2( -1.04, 1.04 ) * a004 + r.zw;\n}\nfloat punctualLightIntensityToIrradianceFactor( const in float lightDistance, const in float cutoffDistance, const in float decayExponent ) {\n#if defined ( PHYSICALLY_CORRECT_LIGHTS )\n\tfloat distanceFalloff = 1.0 / max( pow( lightDistance, decayExponent ), 0.01 );\n\tif( cutoffDistance > 0.0 ) {\n\t\tdistanceFalloff *= pow2( saturate( 1.0 - pow4( lightDistance / cutoffDistance ) ) );\n\t}\n\treturn distanceFalloff;\n#else\n\tif( cutoffDistance > 0.0 && decayExponent > 0.0 ) {\n\t\treturn pow( saturate( -lightDistance / cutoffDistance + 1.0 ), decayExponent );\n\t}\n\treturn 1.0;\n#endif\n}\nvec3 BRDF_Diffuse_Lambert( const in vec3 diffuseColor ) {\n\treturn RECIPROCAL_PI * diffuseColor;\n}\nvec3 F_Schlick( const in vec3 specularColor, const in float dotLH ) {\n\tfloat fresnel = exp2( ( -5.55473 * dotLH - 6.98316 ) * dotLH );\n\treturn ( 1.0 - specularColor ) * fresnel + specularColor;\n}\nvec3 F_Schlick_RoughnessDependent( const in vec3 F0, const in float dotNV, const in float roughness ) {\n\tfloat fresnel = exp2( ( -5.55473 * dotNV - 6.98316 ) * dotNV );\n\tvec3 Fr = max( vec3( 1.0 - roughness ), F0 ) - F0;\n\treturn Fr * fresnel + F0;\n}\nfloat G_GGX_Smith( const in float alpha, const in float dotNL, const in float dotNV ) {\n\tfloat a2 = pow2( alpha );\n\tfloat gl = dotNL + sqrt( a2 + ( 1.0 - a2 ) * pow2( dotNL ) );\n\tfloat gv = dotNV + sqrt( a2 + ( 1.0 - a2 ) * pow2( dotNV ) );\n\treturn 1.0 / ( gl * gv );\n}\nfloat G_GGX_SmithCorrelated( const in float alpha, const in float dotNL, const in float dotNV ) {\n\tfloat a2 = pow2( alpha );\n\tfloat gv = dotNL * sqrt( a2 + ( 1.0 - a2 ) * pow2( dotNV ) );\n\tfloat gl = dotNV * sqrt( a2 + ( 1.0 - a2 ) * pow2( dotNL ) );\n\treturn 0.5 / max( gv + gl, EPSILON );\n}\nfloat D_GGX( const in float alpha, const in float dotNH ) {\n\tfloat a2 = pow2( alpha );\n\tfloat denom = pow2( dotNH ) * ( a2 - 1.0 ) + 1.0;\n\treturn RECIPROCAL_PI * a2 / pow2( denom );\n}\nvec3 BRDF_Specular_GGX( const in IncidentLight incidentLight, const in vec3 viewDir, const in vec3 normal, const in vec3 specularColor, const in float roughness ) {\n\tfloat alpha = pow2( roughness );\n\tvec3 halfDir = normalize( incidentLight.direction + viewDir );\n\tfloat dotNL = saturate( dot( normal, incidentLight.direction ) );\n\tfloat dotNV = saturate( dot( normal, viewDir ) );\n\tfloat dotNH = saturate( dot( normal, halfDir ) );\n\tfloat dotLH = saturate( dot( incidentLight.direction, halfDir ) );\n\tvec3 F = F_Schlick( specularColor, dotLH );\n\tfloat G = G_GGX_SmithCorrelated( alpha, dotNL, dotNV );\n\tfloat D = D_GGX( alpha, dotNH );\n\treturn F * ( G * D );\n}\nvec2 LTC_Uv( const in vec3 N, const in vec3 V, const in float roughness ) {\n\tconst float LUT_SIZE = 64.0;\n\tconst float LUT_SCALE = ( LUT_SIZE - 1.0 ) / LUT_SIZE;\n\tconst float LUT_BIAS = 0.5 / LUT_SIZE;\n\tfloat dotNV = saturate( dot( N, V ) );\n\tvec2 uv = vec2( roughness, sqrt( 1.0 - dotNV ) );\n\tuv = uv * LUT_SCALE + LUT_BIAS;\n\treturn uv;\n}\nfloat LTC_ClippedSphereFormFactor( const in vec3 f ) {\n\tfloat l = length( f );\n\treturn max( ( l * l + f.z ) / ( l + 1.0 ), 0.0 );\n}\nvec3 LTC_EdgeVectorFormFactor( const in vec3 v1, const in vec3 v2 ) {\n\tfloat x = dot( v1, v2 );\n\tfloat y = abs( x );\n\tfloat a = 0.8543985 + ( 0.4965155 + 0.0145206 * y ) * y;\n\tfloat b = 3.4175940 + ( 4.1616724 + y ) * y;\n\tfloat v = a / b;\n\tfloat theta_sintheta = ( x > 0.0 ) ? v : 0.5 * inversesqrt( max( 1.0 - x * x, 1e-7 ) ) - v;\n\treturn cross( v1, v2 ) * theta_sintheta;\n}\nvec3 LTC_Evaluate( const in vec3 N, const in vec3 V, const in vec3 P, const in mat3 mInv, const in vec3 rectCoords[ 4 ] ) {\n\tvec3 v1 = rectCoords[ 1 ] - rectCoords[ 0 ];\n\tvec3 v2 = rectCoords[ 3 ] - rectCoords[ 0 ];\n\tvec3 lightNormal = cross( v1, v2 );\n\tif( dot( lightNormal, P - rectCoords[ 0 ] ) < 0.0 ) return vec3( 0.0 );\n\tvec3 T1, T2;\n\tT1 = normalize( V - N * dot( V, N ) );\n\tT2 = - cross( N, T1 );\n\tmat3 mat = mInv * transposeMat3( mat3( T1, T2, N ) );\n\tvec3 coords[ 4 ];\n\tcoords[ 0 ] = mat * ( rectCoords[ 0 ] - P );\n\tcoords[ 1 ] = mat * ( rectCoords[ 1 ] - P );\n\tcoords[ 2 ] = mat * ( rectCoords[ 2 ] - P );\n\tcoords[ 3 ] = mat * ( rectCoords[ 3 ] - P );\n\tcoords[ 0 ] = normalize( coords[ 0 ] );\n\tcoords[ 1 ] = normalize( coords[ 1 ] );\n\tcoords[ 2 ] = normalize( coords[ 2 ] );\n\tcoords[ 3 ] = normalize( coords[ 3 ] );\n\tvec3 vectorFormFactor = vec3( 0.0 );\n\tvectorFormFactor += LTC_EdgeVectorFormFactor( coords[ 0 ], coords[ 1 ] );\n\tvectorFormFactor += LTC_EdgeVectorFormFactor( coords[ 1 ], coords[ 2 ] );\n\tvectorFormFactor += LTC_EdgeVectorFormFactor( coords[ 2 ], coords[ 3 ] );\n\tvectorFormFactor += LTC_EdgeVectorFormFactor( coords[ 3 ], coords[ 0 ] );\n\tfloat result = LTC_ClippedSphereFormFactor( vectorFormFactor );\n\treturn vec3( result );\n}\nvec3 BRDF_Specular_GGX_Environment( const in vec3 viewDir, const in vec3 normal, const in vec3 specularColor, const in float roughness ) {\n\tfloat dotNV = saturate( dot( normal, viewDir ) );\n\tvec2 brdf = integrateSpecularBRDF( dotNV, roughness );\n\treturn specularColor * brdf.x + brdf.y;\n}\nvoid BRDF_Specular_Multiscattering_Environment( const in GeometricContext geometry, const in vec3 specularColor, const in float roughness, inout vec3 singleScatter, inout vec3 multiScatter ) {\n\tfloat dotNV = saturate( dot( geometry.normal, geometry.viewDir ) );\n\tvec3 F = F_Schlick_RoughnessDependent( specularColor, dotNV, roughness );\n\tvec2 brdf = integrateSpecularBRDF( dotNV, roughness );\n\tvec3 FssEss = F * brdf.x + brdf.y;\n\tfloat Ess = brdf.x + brdf.y;\n\tfloat Ems = 1.0 - Ess;\n\tvec3 Favg = specularColor + ( 1.0 - specularColor ) * 0.047619;\tvec3 Fms = FssEss * Favg / ( 1.0 - Ems * Favg );\n\tsingleScatter += FssEss;\n\tmultiScatter += Fms * Ems;\n}\nfloat G_BlinnPhong_Implicit( ) {\n\treturn 0.25;\n}\nfloat D_BlinnPhong( const in float shininess, const in float dotNH ) {\n\treturn RECIPROCAL_PI * ( shininess * 0.5 + 1.0 ) * pow( dotNH, shininess );\n}\nvec3 BRDF_Specular_BlinnPhong( const in IncidentLight incidentLight, const in GeometricContext geometry, const in vec3 specularColor, const in float shininess ) {\n\tvec3 halfDir = normalize( incidentLight.direction + geometry.viewDir );\n\tfloat dotNH = saturate( dot( geometry.normal, halfDir ) );\n\tfloat dotLH = saturate( dot( incidentLight.direction, halfDir ) );\n\tvec3 F = F_Schlick( specularColor, dotLH );\n\tfloat G = G_BlinnPhong_Implicit( );\n\tfloat D = D_BlinnPhong( shininess, dotNH );\n\treturn F * ( G * D );\n}\nfloat GGXRoughnessToBlinnExponent( const in float ggxRoughness ) {\n\treturn ( 2.0 / pow2( ggxRoughness + 0.0001 ) - 2.0 );\n}\nfloat BlinnExponentToGGXRoughness( const in float blinnExponent ) {\n\treturn sqrt( 2.0 / ( blinnExponent + 2.0 ) );\n}\n#if defined( USE_SHEEN )\nfloat D_Charlie(float roughness, float NoH) {\n\tfloat invAlpha = 1.0 / roughness;\n\tfloat cos2h = NoH * NoH;\n\tfloat sin2h = max(1.0 - cos2h, 0.0078125);\treturn (2.0 + invAlpha) * pow(sin2h, invAlpha * 0.5) / (2.0 * PI);\n}\nfloat V_Neubelt(float NoV, float NoL) {\n\treturn saturate(1.0 / (4.0 * (NoL + NoV - NoL * NoV)));\n}\nvec3 BRDF_Specular_Sheen( const in float roughness, const in vec3 L, const in GeometricContext geometry, vec3 specularColor ) {\n\tvec3 N = geometry.normal;\n\tvec3 V = geometry.viewDir;\n\tvec3 H = normalize( V + L );\n\tfloat dotNH = saturate( dot( N, H ) );\n\treturn specularColor * D_Charlie( roughness, dotNH ) * V_Neubelt( dot(N, V), dot(N, L) );\n}\n#endif",
	bumpmap_pars_fragment = "#ifdef USE_BUMPMAP\n\tuniform sampler2D bumpMap;\n\tuniform float bumpScale;\n\tvec2 dHdxy_fwd() {\n\t\tvec2 dSTdx = dFdx( vUv );\n\t\tvec2 dSTdy = dFdy( vUv );\n\t\tfloat Hll = bumpScale * texture2D( bumpMap, vUv ).x;\n\t\tfloat dBx = bumpScale * texture2D( bumpMap, vUv + dSTdx ).x - Hll;\n\t\tfloat dBy = bumpScale * texture2D( bumpMap, vUv + dSTdy ).x - Hll;\n\t\treturn vec2( dBx, dBy );\n\t}\n\tvec3 perturbNormalArb( vec3 surf_pos, vec3 surf_norm, vec2 dHdxy, float faceDirection ) {\n\t\tvec3 vSigmaX = vec3( dFdx( surf_pos.x ), dFdx( surf_pos.y ), dFdx( surf_pos.z ) );\n\t\tvec3 vSigmaY = vec3( dFdy( surf_pos.x ), dFdy( surf_pos.y ), dFdy( surf_pos.z ) );\n\t\tvec3 vN = surf_norm;\n\t\tvec3 R1 = cross( vSigmaY, vN );\n\t\tvec3 R2 = cross( vN, vSigmaX );\n\t\tfloat fDet = dot( vSigmaX, R1 ) * faceDirection;\n\t\tvec3 vGrad = sign( fDet ) * ( dHdxy.x * R1 + dHdxy.y * R2 );\n\t\treturn normalize( abs( fDet ) * surf_norm - vGrad );\n\t}\n#endif",
	clipping_planes_fragment = "#if NUM_CLIPPING_PLANES > 0\n\tvec4 plane;\n\t#pragma unroll_loop_start\n\tfor ( int i = 0; i < UNION_CLIPPING_PLANES; i ++ ) {\n\t\tplane = clippingPlanes[ i ];\n\t\tif ( dot( vClipPosition, plane.xyz ) > plane.w ) discard;\n\t}\n\t#pragma unroll_loop_end\n\t#if UNION_CLIPPING_PLANES < NUM_CLIPPING_PLANES\n\t\tbool clipped = true;\n\t\t#pragma unroll_loop_start\n\t\tfor ( int i = UNION_CLIPPING_PLANES; i < NUM_CLIPPING_PLANES; i ++ ) {\n\t\t\tplane = clippingPlanes[ i ];\n\t\t\tclipped = ( dot( vClipPosition, plane.xyz ) > plane.w ) && clipped;\n\t\t}\n\t\t#pragma unroll_loop_end\n\t\tif ( clipped ) discard;\n\t#endif\n#endif",
	clipping_planes_pars_fragment = "#if NUM_CLIPPING_PLANES > 0\n\tvarying vec3 vClipPosition;\n\tuniform vec4 clippingPlanes[ NUM_CLIPPING_PLANES ];\n#endif",
	clipping_planes_pars_vertex = "#if NUM_CLIPPING_PLANES > 0\n\tvarying vec3 vClipPosition;\n#endif",
	clipping_planes_vertex = "#if NUM_CLIPPING_PLANES > 0\n\tvClipPosition = - mvPosition.xyz;\n#endif",
	color_fragment = "#ifdef USE_COLOR\n\tdiffuseColor.rgb *= vColor;\n#endif",
	color_pars_fragment = "#ifdef USE_COLOR\n\tvarying vec3 vColor;\n#endif",
	color_pars_vertex = "#if defined( USE_COLOR ) || defined( USE_INSTANCING_COLOR )\n\tvarying vec3 vColor;\n#endif",
	color_vertex = "#if defined( USE_COLOR ) || defined( USE_INSTANCING_COLOR )\n\tvColor = vec3( 1.0 );\n#endif\n#ifdef USE_COLOR\n\tvColor.xyz *= color.xyz;\n#endif\n#ifdef USE_INSTANCING_COLOR\n\tvColor.xyz *= instanceColor.xyz;\n#endif",
	common = "#define PI 3.141592653589793\n#define PI2 6.283185307179586\n#define PI_HALF 1.5707963267948966\n#define RECIPROCAL_PI 0.3183098861837907\n#define RECIPROCAL_PI2 0.15915494309189535\n#define EPSILON 1e-6\n#ifndef saturate\n#define saturate(a) clamp( a, 0.0, 1.0 )\n#endif\n#define whiteComplement(a) ( 1.0 - saturate( a ) )\nfloat pow2( const in float x ) { return x*x; }\nfloat pow3( const in float x ) { return x*x*x; }\nfloat pow4( const in float x ) { float x2 = x*x; return x2*x2; }\nfloat average( const in vec3 color ) { return dot( color, vec3( 0.3333 ) ); }\nhighp float rand( const in vec2 uv ) {\n\tconst highp float a = 12.9898, b = 78.233, c = 43758.5453;\n\thighp float dt = dot( uv.xy, vec2( a,b ) ), sn = mod( dt, PI );\n\treturn fract(sin(sn) * c);\n}\n#ifdef HIGH_PRECISION\n\tfloat precisionSafeLength( vec3 v ) { return length( v ); }\n#else\n\tfloat max3( vec3 v ) { return max( max( v.x, v.y ), v.z ); }\n\tfloat precisionSafeLength( vec3 v ) {\n\t\tfloat maxComponent = max3( abs( v ) );\n\t\treturn length( v / maxComponent ) * maxComponent;\n\t}\n#endif\nstruct IncidentLight {\n\tvec3 color;\n\tvec3 direction;\n\tbool visible;\n};\nstruct ReflectedLight {\n\tvec3 directDiffuse;\n\tvec3 directSpecular;\n\tvec3 indirectDiffuse;\n\tvec3 indirectSpecular;\n};\nstruct GeometricContext {\n\tvec3 position;\n\tvec3 normal;\n\tvec3 viewDir;\n#ifdef CLEARCOAT\n\tvec3 clearcoatNormal;\n#endif\n};\nvec3 transformDirection( in vec3 dir, in mat4 matrix ) {\n\treturn normalize( ( matrix * vec4( dir, 0.0 ) ).xyz );\n}\nvec3 inverseTransformDirection( in vec3 dir, in mat4 matrix ) {\n\treturn normalize( ( vec4( dir, 0.0 ) * matrix ).xyz );\n}\nvec3 projectOnPlane(in vec3 point, in vec3 pointOnPlane, in vec3 planeNormal ) {\n\tfloat distance = dot( planeNormal, point - pointOnPlane );\n\treturn - distance * planeNormal + point;\n}\nfloat sideOfPlane( in vec3 point, in vec3 pointOnPlane, in vec3 planeNormal ) {\n\treturn sign( dot( point - pointOnPlane, planeNormal ) );\n}\nvec3 linePlaneIntersect( in vec3 pointOnLine, in vec3 lineDirection, in vec3 pointOnPlane, in vec3 planeNormal ) {\n\treturn lineDirection * ( dot( planeNormal, pointOnPlane - pointOnLine ) / dot( planeNormal, lineDirection ) ) + pointOnLine;\n}\nmat3 transposeMat3( const in mat3 m ) {\n\tmat3 tmp;\n\ttmp[ 0 ] = vec3( m[ 0 ].x, m[ 1 ].x, m[ 2 ].x );\n\ttmp[ 1 ] = vec3( m[ 0 ].y, m[ 1 ].y, m[ 2 ].y );\n\ttmp[ 2 ] = vec3( m[ 0 ].z, m[ 1 ].z, m[ 2 ].z );\n\treturn tmp;\n}\nfloat linearToRelativeLuminance( const in vec3 color ) {\n\tvec3 weights = vec3( 0.2126, 0.7152, 0.0722 );\n\treturn dot( weights, color.rgb );\n}\nbool isPerspectiveMatrix( mat4 m ) {\n\treturn m[ 2 ][ 3 ] == - 1.0;\n}\nvec2 equirectUv( in vec3 dir ) {\n\tfloat u = atan( dir.z, dir.x ) * RECIPROCAL_PI2 + 0.5;\n\tfloat v = asin( clamp( dir.y, - 1.0, 1.0 ) ) * RECIPROCAL_PI + 0.5;\n\treturn vec2( u, v );\n}",
	cube_uv_reflection_fragment = "#ifdef ENVMAP_TYPE_CUBE_UV\n\t#define cubeUV_maxMipLevel 8.0\n\t#define cubeUV_minMipLevel 4.0\n\t#define cubeUV_maxTileSize 256.0\n\t#define cubeUV_minTileSize 16.0\n\tfloat getFace( vec3 direction ) {\n\t\tvec3 absDirection = abs( direction );\n\t\tfloat face = - 1.0;\n\t\tif ( absDirection.x > absDirection.z ) {\n\t\t\tif ( absDirection.x > absDirection.y )\n\t\t\t\tface = direction.x > 0.0 ? 0.0 : 3.0;\n\t\t\telse\n\t\t\t\tface = direction.y > 0.0 ? 1.0 : 4.0;\n\t\t} else {\n\t\t\tif ( absDirection.z > absDirection.y )\n\t\t\t\tface = direction.z > 0.0 ? 2.0 : 5.0;\n\t\t\telse\n\t\t\t\tface = direction.y > 0.0 ? 1.0 : 4.0;\n\t\t}\n\t\treturn face;\n\t}\n\tvec2 getUV( vec3 direction, float face ) {\n\t\tvec2 uv;\n\t\tif ( face == 0.0 ) {\n\t\t\tuv = vec2( direction.z, direction.y ) / abs( direction.x );\n\t\t} else if ( face == 1.0 ) {\n\t\t\tuv = vec2( - direction.x, - direction.z ) / abs( direction.y );\n\t\t} else if ( face == 2.0 ) {\n\t\t\tuv = vec2( - direction.x, direction.y ) / abs( direction.z );\n\t\t} else if ( face == 3.0 ) {\n\t\t\tuv = vec2( - direction.z, direction.y ) / abs( direction.x );\n\t\t} else if ( face == 4.0 ) {\n\t\t\tuv = vec2( - direction.x, direction.z ) / abs( direction.y );\n\t\t} else {\n\t\t\tuv = vec2( direction.x, direction.y ) / abs( direction.z );\n\t\t}\n\t\treturn 0.5 * ( uv + 1.0 );\n\t}\n\tvec3 bilinearCubeUV( sampler2D envMap, vec3 direction, float mipInt ) {\n\t\tfloat face = getFace( direction );\n\t\tfloat filterInt = max( cubeUV_minMipLevel - mipInt, 0.0 );\n\t\tmipInt = max( mipInt, cubeUV_minMipLevel );\n\t\tfloat faceSize = exp2( mipInt );\n\t\tfloat texelSize = 1.0 / ( 3.0 * cubeUV_maxTileSize );\n\t\tvec2 uv = getUV( direction, face ) * ( faceSize - 1.0 );\n\t\tvec2 f = fract( uv );\n\t\tuv += 0.5 - f;\n\t\tif ( face > 2.0 ) {\n\t\t\tuv.y += faceSize;\n\t\t\tface -= 3.0;\n\t\t}\n\t\tuv.x += face * faceSize;\n\t\tif ( mipInt < cubeUV_maxMipLevel ) {\n\t\t\tuv.y += 2.0 * cubeUV_maxTileSize;\n\t\t}\n\t\tuv.y += filterInt * 2.0 * cubeUV_minTileSize;\n\t\tuv.x += 3.0 * max( 0.0, cubeUV_maxTileSize - 2.0 * faceSize );\n\t\tuv *= texelSize;\n\t\tvec3 tl = envMapTexelToLinear( texture2D( envMap, uv ) ).rgb;\n\t\tuv.x += texelSize;\n\t\tvec3 tr = envMapTexelToLinear( texture2D( envMap, uv ) ).rgb;\n\t\tuv.y += texelSize;\n\t\tvec3 br = envMapTexelToLinear( texture2D( envMap, uv ) ).rgb;\n\t\tuv.x -= texelSize;\n\t\tvec3 bl = envMapTexelToLinear( texture2D( envMap, uv ) ).rgb;\n\t\tvec3 tm = mix( tl, tr, f.x );\n\t\tvec3 bm = mix( bl, br, f.x );\n\t\treturn mix( tm, bm, f.y );\n\t}\n\t#define r0 1.0\n\t#define v0 0.339\n\t#define m0 - 2.0\n\t#define r1 0.8\n\t#define v1 0.276\n\t#define m1 - 1.0\n\t#define r4 0.4\n\t#define v4 0.046\n\t#define m4 2.0\n\t#define r5 0.305\n\t#define v5 0.016\n\t#define m5 3.0\n\t#define r6 0.21\n\t#define v6 0.0038\n\t#define m6 4.0\n\tfloat roughnessToMip( float roughness ) {\n\t\tfloat mip = 0.0;\n\t\tif ( roughness >= r1 ) {\n\t\t\tmip = ( r0 - roughness ) * ( m1 - m0 ) / ( r0 - r1 ) + m0;\n\t\t} else if ( roughness >= r4 ) {\n\t\t\tmip = ( r1 - roughness ) * ( m4 - m1 ) / ( r1 - r4 ) + m1;\n\t\t} else if ( roughness >= r5 ) {\n\t\t\tmip = ( r4 - roughness ) * ( m5 - m4 ) / ( r4 - r5 ) + m4;\n\t\t} else if ( roughness >= r6 ) {\n\t\t\tmip = ( r5 - roughness ) * ( m6 - m5 ) / ( r5 - r6 ) + m5;\n\t\t} else {\n\t\t\tmip = - 2.0 * log2( 1.16 * roughness );\t\t}\n\t\treturn mip;\n\t}\n\tvec4 textureCubeUV( sampler2D envMap, vec3 sampleDir, float roughness ) {\n\t\tfloat mip = clamp( roughnessToMip( roughness ), m0, cubeUV_maxMipLevel );\n\t\tfloat mipF = fract( mip );\n\t\tfloat mipInt = floor( mip );\n\t\tvec3 color0 = bilinearCubeUV( envMap, sampleDir, mipInt );\n\t\tif ( mipF == 0.0 ) {\n\t\t\treturn vec4( color0, 1.0 );\n\t\t} else {\n\t\t\tvec3 color1 = bilinearCubeUV( envMap, sampleDir, mipInt + 1.0 );\n\t\t\treturn vec4( mix( color0, color1, mipF ), 1.0 );\n\t\t}\n\t}\n#endif",
	defaultnormal_vertex = "vec3 transformedNormal = objectNormal;\n#ifdef USE_INSTANCING\n\tmat3 m = mat3( instanceMatrix );\n\ttransformedNormal /= vec3( dot( m[ 0 ], m[ 0 ] ), dot( m[ 1 ], m[ 1 ] ), dot( m[ 2 ], m[ 2 ] ) );\n\ttransformedNormal = m * transformedNormal;\n#endif\ntransformedNormal = normalMatrix * transformedNormal;\n#ifdef FLIP_SIDED\n\ttransformedNormal = - transformedNormal;\n#endif\n#ifdef USE_TANGENT\n\tvec3 transformedTangent = ( modelViewMatrix * vec4( objectTangent, 0.0 ) ).xyz;\n\t#ifdef FLIP_SIDED\n\t\ttransformedTangent = - transformedTangent;\n\t#endif\n#endif",
	displacementmap_pars_vertex = "#ifdef USE_DISPLACEMENTMAP\n\tuniform sampler2D displacementMap;\n\tuniform float displacementScale;\n\tuniform float displacementBias;\n#endif",
	displacementmap_vertex = "#ifdef USE_DISPLACEMENTMAP\n\ttransformed += normalize( objectNormal ) * ( texture2D( displacementMap, vUv ).x * displacementScale + displacementBias );\n#endif",
	emissivemap_fragment = "#ifdef USE_EMISSIVEMAP\n\tvec4 emissiveColor = texture2D( emissiveMap, vUv );\n\temissiveColor.rgb = emissiveMapTexelToLinear( emissiveColor ).rgb;\n\ttotalEmissiveRadiance *= emissiveColor.rgb;\n#endif",
	emissivemap_pars_fragment = "#ifdef USE_EMISSIVEMAP\n\tuniform sampler2D emissiveMap;\n#endif",
	encodings_fragment = "gl_FragColor = linearToOutputTexel( gl_FragColor );",
	encodings_pars_fragment = "\nvec4 LinearToLinear( in vec4 value ) {\n\treturn value;\n}\nvec4 GammaToLinear( in vec4 value, in float gammaFactor ) {\n\treturn vec4( pow( value.rgb, vec3( gammaFactor ) ), value.a );\n}\nvec4 LinearToGamma( in vec4 value, in float gammaFactor ) {\n\treturn vec4( pow( value.rgb, vec3( 1.0 / gammaFactor ) ), value.a );\n}\nvec4 sRGBToLinear( in vec4 value ) {\n\treturn vec4( mix( pow( value.rgb * 0.9478672986 + vec3( 0.0521327014 ), vec3( 2.4 ) ), value.rgb * 0.0773993808, vec3( lessThanEqual( value.rgb, vec3( 0.04045 ) ) ) ), value.a );\n}\nvec4 LinearTosRGB( in vec4 value ) {\n\treturn vec4( mix( pow( value.rgb, vec3( 0.41666 ) ) * 1.055 - vec3( 0.055 ), value.rgb * 12.92, vec3( lessThanEqual( value.rgb, vec3( 0.0031308 ) ) ) ), value.a );\n}\nvec4 RGBEToLinear( in vec4 value ) {\n\treturn vec4( value.rgb * exp2( value.a * 255.0 - 128.0 ), 1.0 );\n}\nvec4 LinearToRGBE( in vec4 value ) {\n\tfloat maxComponent = max( max( value.r, value.g ), value.b );\n\tfloat fExp = clamp( ceil( log2( maxComponent ) ), -128.0, 127.0 );\n\treturn vec4( value.rgb / exp2( fExp ), ( fExp + 128.0 ) / 255.0 );\n}\nvec4 RGBMToLinear( in vec4 value, in float maxRange ) {\n\treturn vec4( value.rgb * value.a * maxRange, 1.0 );\n}\nvec4 LinearToRGBM( in vec4 value, in float maxRange ) {\n\tfloat maxRGB = max( value.r, max( value.g, value.b ) );\n\tfloat M = clamp( maxRGB / maxRange, 0.0, 1.0 );\n\tM = ceil( M * 255.0 ) / 255.0;\n\treturn vec4( value.rgb / ( M * maxRange ), M );\n}\nvec4 RGBDToLinear( in vec4 value, in float maxRange ) {\n\treturn vec4( value.rgb * ( ( maxRange / 255.0 ) / value.a ), 1.0 );\n}\nvec4 LinearToRGBD( in vec4 value, in float maxRange ) {\n\tfloat maxRGB = max( value.r, max( value.g, value.b ) );\n\tfloat D = max( maxRange / maxRGB, 1.0 );\n\tD = clamp( floor( D ) / 255.0, 0.0, 1.0 );\n\treturn vec4( value.rgb * ( D * ( 255.0 / maxRange ) ), D );\n}\nconst mat3 cLogLuvM = mat3( 0.2209, 0.3390, 0.4184, 0.1138, 0.6780, 0.7319, 0.0102, 0.1130, 0.2969 );\nvec4 LinearToLogLuv( in vec4 value ) {\n\tvec3 Xp_Y_XYZp = cLogLuvM * value.rgb;\n\tXp_Y_XYZp = max( Xp_Y_XYZp, vec3( 1e-6, 1e-6, 1e-6 ) );\n\tvec4 vResult;\n\tvResult.xy = Xp_Y_XYZp.xy / Xp_Y_XYZp.z;\n\tfloat Le = 2.0 * log2(Xp_Y_XYZp.y) + 127.0;\n\tvResult.w = fract( Le );\n\tvResult.z = ( Le - ( floor( vResult.w * 255.0 ) ) / 255.0 ) / 255.0;\n\treturn vResult;\n}\nconst mat3 cLogLuvInverseM = mat3( 6.0014, -2.7008, -1.7996, -1.3320, 3.1029, -5.7721, 0.3008, -1.0882, 5.6268 );\nvec4 LogLuvToLinear( in vec4 value ) {\n\tfloat Le = value.z * 255.0 + value.w;\n\tvec3 Xp_Y_XYZp;\n\tXp_Y_XYZp.y = exp2( ( Le - 127.0 ) / 2.0 );\n\tXp_Y_XYZp.z = Xp_Y_XYZp.y / value.y;\n\tXp_Y_XYZp.x = value.x * Xp_Y_XYZp.z;\n\tvec3 vRGB = cLogLuvInverseM * Xp_Y_XYZp.rgb;\n\treturn vec4( max( vRGB, 0.0 ), 1.0 );\n}",
	envmap_fragment = "#ifdef USE_ENVMAP\n\t#ifdef ENV_WORLDPOS\n\t\tvec3 cameraToFrag;\n\t\tif ( isOrthographic ) {\n\t\t\tcameraToFrag = normalize( vec3( - viewMatrix[ 0 ][ 2 ], - viewMatrix[ 1 ][ 2 ], - viewMatrix[ 2 ][ 2 ] ) );\n\t\t} else {\n\t\t\tcameraToFrag = normalize( vWorldPosition - cameraPosition );\n\t\t}\n\t\tvec3 worldNormal = inverseTransformDirection( normal, viewMatrix );\n\t\t#ifdef ENVMAP_MODE_REFLECTION\n\t\t\tvec3 reflectVec = reflect( cameraToFrag, worldNormal );\n\t\t#else\n\t\t\tvec3 reflectVec = refract( cameraToFrag, worldNormal, refractionRatio );\n\t\t#endif\n\t#else\n\t\tvec3 reflectVec = vReflect;\n\t#endif\n\t#ifdef ENVMAP_TYPE_CUBE\n\t\tvec4 envColor = textureCube( envMap, vec3( flipEnvMap * reflectVec.x, reflectVec.yz ) );\n\t#elif defined( ENVMAP_TYPE_CUBE_UV )\n\t\tvec4 envColor = textureCubeUV( envMap, reflectVec, 0.0 );\n\t#else\n\t\tvec4 envColor = vec4( 0.0 );\n\t#endif\n\t#ifndef ENVMAP_TYPE_CUBE_UV\n\t\tenvColor = envMapTexelToLinear( envColor );\n\t#endif\n\t#ifdef ENVMAP_BLENDING_MULTIPLY\n\t\toutgoingLight = mix( outgoingLight, outgoingLight * envColor.xyz, specularStrength * reflectivity );\n\t#elif defined( ENVMAP_BLENDING_MIX )\n\t\toutgoingLight = mix( outgoingLight, envColor.xyz, specularStrength * reflectivity );\n\t#elif defined( ENVMAP_BLENDING_ADD )\n\t\toutgoingLight += envColor.xyz * specularStrength * reflectivity;\n\t#endif\n#endif",
	envmap_common_pars_fragment = "#ifdef USE_ENVMAP\n\tuniform float envMapIntensity;\n\tuniform float flipEnvMap;\n\tuniform int maxMipLevel;\n\t#ifdef ENVMAP_TYPE_CUBE\n\t\tuniform samplerCube envMap;\n\t#else\n\t\tuniform sampler2D envMap;\n\t#endif\n\t\n#endif",
	envmap_pars_fragment = "#ifdef USE_ENVMAP\n\tuniform float reflectivity;\n\t#if defined( USE_BUMPMAP ) || defined( USE_NORMALMAP ) || defined( PHONG )\n\t\t#define ENV_WORLDPOS\n\t#endif\n\t#ifdef ENV_WORLDPOS\n\t\tvarying vec3 vWorldPosition;\n\t\tuniform float refractionRatio;\n\t#else\n\t\tvarying vec3 vReflect;\n\t#endif\n#endif",
	envmap_pars_vertex = "#ifdef USE_ENVMAP\n\t#if defined( USE_BUMPMAP ) || defined( USE_NORMALMAP ) ||defined( PHONG )\n\t\t#define ENV_WORLDPOS\n\t#endif\n\t#ifdef ENV_WORLDPOS\n\t\t\n\t\tvarying vec3 vWorldPosition;\n\t#else\n\t\tvarying vec3 vReflect;\n\t\tuniform float refractionRatio;\n\t#endif\n#endif",
	envmap_vertex = "#ifdef USE_ENVMAP\n\t#ifdef ENV_WORLDPOS\n\t\tvWorldPosition = worldPosition.xyz;\n\t#else\n\t\tvec3 cameraToVertex;\n\t\tif ( isOrthographic ) {\n\t\t\tcameraToVertex = normalize( vec3( - viewMatrix[ 0 ][ 2 ], - viewMatrix[ 1 ][ 2 ], - viewMatrix[ 2 ][ 2 ] ) );\n\t\t} else {\n\t\t\tcameraToVertex = normalize( worldPosition.xyz - cameraPosition );\n\t\t}\n\t\tvec3 worldNormal = inverseTransformDirection( transformedNormal, viewMatrix );\n\t\t#ifdef ENVMAP_MODE_REFLECTION\n\t\t\tvReflect = reflect( cameraToVertex, worldNormal );\n\t\t#else\n\t\t\tvReflect = refract( cameraToVertex, worldNormal, refractionRatio );\n\t\t#endif\n\t#endif\n#endif",
	fog_vertex = "#ifdef USE_FOG\n\tfogDepth = - mvPosition.z;\n#endif",
	fog_pars_vertex = "#ifdef USE_FOG\n\tvarying float fogDepth;\n#endif",
	fog_fragment = "#ifdef USE_FOG\n\t#ifdef FOG_EXP2\n\t\tfloat fogFactor = 1.0 - exp( - fogDensity * fogDensity * fogDepth * fogDepth );\n\t#else\n\t\tfloat fogFactor = smoothstep( fogNear, fogFar, fogDepth );\n\t#endif\n\tgl_FragColor.rgb = mix( gl_FragColor.rgb, fogColor, fogFactor );\n#endif",
	fog_pars_fragment = "#ifdef USE_FOG\n\tuniform vec3 fogColor;\n\tvarying float fogDepth;\n\t#ifdef FOG_EXP2\n\t\tuniform float fogDensity;\n\t#else\n\t\tuniform float fogNear;\n\t\tuniform float fogFar;\n\t#endif\n#endif",
	gradientmap_pars_fragment = "#ifdef USE_GRADIENTMAP\n\tuniform sampler2D gradientMap;\n#endif\nvec3 getGradientIrradiance( vec3 normal, vec3 lightDirection ) {\n\tfloat dotNL = dot( normal, lightDirection );\n\tvec2 coord = vec2( dotNL * 0.5 + 0.5, 0.0 );\n\t#ifdef USE_GRADIENTMAP\n\t\treturn texture2D( gradientMap, coord ).rgb;\n\t#else\n\t\treturn ( coord.x < 0.7 ) ? vec3( 0.7 ) : vec3( 1.0 );\n\t#endif\n}",
	lightmap_fragment = "#ifdef USE_LIGHTMAP\n\tvec4 lightMapTexel= texture2D( lightMap, vUv2 );\n\treflectedLight.indirectDiffuse += PI * lightMapTexelToLinear( lightMapTexel ).rgb * lightMapIntensity;\n#endif",
	lightmap_pars_fragment = "#ifdef USE_LIGHTMAP\n\tuniform sampler2D lightMap;\n\tuniform float lightMapIntensity;\n#endif",
	lights_lambert_vertex = "vec3 diffuse = vec3( 1.0 );\nGeometricContext geometry;\ngeometry.position = mvPosition.xyz;\ngeometry.normal = normalize( transformedNormal );\ngeometry.viewDir = ( isOrthographic ) ? vec3( 0, 0, 1 ) : normalize( -mvPosition.xyz );\nGeometricContext backGeometry;\nbackGeometry.position = geometry.position;\nbackGeometry.normal = -geometry.normal;\nbackGeometry.viewDir = geometry.viewDir;\nvLightFront = vec3( 0.0 );\nvIndirectFront = vec3( 0.0 );\n#ifdef DOUBLE_SIDED\n\tvLightBack = vec3( 0.0 );\n\tvIndirectBack = vec3( 0.0 );\n#endif\nIncidentLight directLight;\nfloat dotNL;\nvec3 directLightColor_Diffuse;\nvIndirectFront += getAmbientLightIrradiance( ambientLightColor );\nvIndirectFront += getLightProbeIrradiance( lightProbe, geometry );\n#ifdef DOUBLE_SIDED\n\tvIndirectBack += getAmbientLightIrradiance( ambientLightColor );\n\tvIndirectBack += getLightProbeIrradiance( lightProbe, backGeometry );\n#endif\n#if NUM_POINT_LIGHTS > 0\n\t#pragma unroll_loop_start\n\tfor ( int i = 0; i < NUM_POINT_LIGHTS; i ++ ) {\n\t\tgetPointDirectLightIrradiance( pointLights[ i ], geometry, directLight );\n\t\tdotNL = dot( geometry.normal, directLight.direction );\n\t\tdirectLightColor_Diffuse = PI * directLight.color;\n\t\tvLightFront += saturate( dotNL ) * directLightColor_Diffuse;\n\t\t#ifdef DOUBLE_SIDED\n\t\t\tvLightBack += saturate( -dotNL ) * directLightColor_Diffuse;\n\t\t#endif\n\t}\n\t#pragma unroll_loop_end\n#endif\n#if NUM_SPOT_LIGHTS > 0\n\t#pragma unroll_loop_start\n\tfor ( int i = 0; i < NUM_SPOT_LIGHTS; i ++ ) {\n\t\tgetSpotDirectLightIrradiance( spotLights[ i ], geometry, directLight );\n\t\tdotNL = dot( geometry.normal, directLight.direction );\n\t\tdirectLightColor_Diffuse = PI * directLight.color;\n\t\tvLightFront += saturate( dotNL ) * directLightColor_Diffuse;\n\t\t#ifdef DOUBLE_SIDED\n\t\t\tvLightBack += saturate( -dotNL ) * directLightColor_Diffuse;\n\t\t#endif\n\t}\n\t#pragma unroll_loop_end\n#endif\n#if NUM_DIR_LIGHTS > 0\n\t#pragma unroll_loop_start\n\tfor ( int i = 0; i < NUM_DIR_LIGHTS; i ++ ) {\n\t\tgetDirectionalDirectLightIrradiance( directionalLights[ i ], geometry, directLight );\n\t\tdotNL = dot( geometry.normal, directLight.direction );\n\t\tdirectLightColor_Diffuse = PI * directLight.color;\n\t\tvLightFront += saturate( dotNL ) * directLightColor_Diffuse;\n\t\t#ifdef DOUBLE_SIDED\n\t\t\tvLightBack += saturate( -dotNL ) * directLightColor_Diffuse;\n\t\t#endif\n\t}\n\t#pragma unroll_loop_end\n#endif\n#if NUM_HEMI_LIGHTS > 0\n\t#pragma unroll_loop_start\n\tfor ( int i = 0; i < NUM_HEMI_LIGHTS; i ++ ) {\n\t\tvIndirectFront += getHemisphereLightIrradiance( hemisphereLights[ i ], geometry );\n\t\t#ifdef DOUBLE_SIDED\n\t\t\tvIndirectBack += getHemisphereLightIrradiance( hemisphereLights[ i ], backGeometry );\n\t\t#endif\n\t}\n\t#pragma unroll_loop_end\n#endif",
	lights_pars_begin = "uniform bool receiveShadow;\nuniform vec3 ambientLightColor;\nuniform vec3 lightProbe[ 9 ];\nvec3 shGetIrradianceAt( in vec3 normal, in vec3 shCoefficients[ 9 ] ) {\n\tfloat x = normal.x, y = normal.y, z = normal.z;\n\tvec3 result = shCoefficients[ 0 ] * 0.886227;\n\tresult += shCoefficients[ 1 ] * 2.0 * 0.511664 * y;\n\tresult += shCoefficients[ 2 ] * 2.0 * 0.511664 * z;\n\tresult += shCoefficients[ 3 ] * 2.0 * 0.511664 * x;\n\tresult += shCoefficients[ 4 ] * 2.0 * 0.429043 * x * y;\n\tresult += shCoefficients[ 5 ] * 2.0 * 0.429043 * y * z;\n\tresult += shCoefficients[ 6 ] * ( 0.743125 * z * z - 0.247708 );\n\tresult += shCoefficients[ 7 ] * 2.0 * 0.429043 * x * z;\n\tresult += shCoefficients[ 8 ] * 0.429043 * ( x * x - y * y );\n\treturn result;\n}\nvec3 getLightProbeIrradiance( const in vec3 lightProbe[ 9 ], const in GeometricContext geometry ) {\n\tvec3 worldNormal = inverseTransformDirection( geometry.normal, viewMatrix );\n\tvec3 irradiance = shGetIrradianceAt( worldNormal, lightProbe );\n\treturn irradiance;\n}\nvec3 getAmbientLightIrradiance( const in vec3 ambientLightColor ) {\n\tvec3 irradiance = ambientLightColor;\n\t#ifndef PHYSICALLY_CORRECT_LIGHTS\n\t\tirradiance *= PI;\n\t#endif\n\treturn irradiance;\n}\n#if NUM_DIR_LIGHTS > 0\n\tstruct DirectionalLight {\n\t\tvec3 direction;\n\t\tvec3 color;\n\t};\n\tuniform DirectionalLight directionalLights[ NUM_DIR_LIGHTS ];\n\tvoid getDirectionalDirectLightIrradiance( const in DirectionalLight directionalLight, const in GeometricContext geometry, out IncidentLight directLight ) {\n\t\tdirectLight.color = directionalLight.color;\n\t\tdirectLight.direction = directionalLight.direction;\n\t\tdirectLight.visible = true;\n\t}\n#endif\n#if NUM_POINT_LIGHTS > 0\n\tstruct PointLight {\n\t\tvec3 position;\n\t\tvec3 color;\n\t\tfloat distance;\n\t\tfloat decay;\n\t};\n\tuniform PointLight pointLights[ NUM_POINT_LIGHTS ];\n\tvoid getPointDirectLightIrradiance( const in PointLight pointLight, const in GeometricContext geometry, out IncidentLight directLight ) {\n\t\tvec3 lVector = pointLight.position - geometry.position;\n\t\tdirectLight.direction = normalize( lVector );\n\t\tfloat lightDistance = length( lVector );\n\t\tdirectLight.color = pointLight.color;\n\t\tdirectLight.color *= punctualLightIntensityToIrradianceFactor( lightDistance, pointLight.distance, pointLight.decay );\n\t\tdirectLight.visible = ( directLight.color != vec3( 0.0 ) );\n\t}\n#endif\n#if NUM_SPOT_LIGHTS > 0\n\tstruct SpotLight {\n\t\tvec3 position;\n\t\tvec3 direction;\n\t\tvec3 color;\n\t\tfloat distance;\n\t\tfloat decay;\n\t\tfloat coneCos;\n\t\tfloat penumbraCos;\n\t};\n\tuniform SpotLight spotLights[ NUM_SPOT_LIGHTS ];\n\tvoid getSpotDirectLightIrradiance( const in SpotLight spotLight, const in GeometricContext geometry, out IncidentLight directLight ) {\n\t\tvec3 lVector = spotLight.position - geometry.position;\n\t\tdirectLight.direction = normalize( lVector );\n\t\tfloat lightDistance = length( lVector );\n\t\tfloat angleCos = dot( directLight.direction, spotLight.direction );\n\t\tif ( angleCos > spotLight.coneCos ) {\n\t\t\tfloat spotEffect = smoothstep( spotLight.coneCos, spotLight.penumbraCos, angleCos );\n\t\t\tdirectLight.color = spotLight.color;\n\t\t\tdirectLight.color *= spotEffect * punctualLightIntensityToIrradianceFactor( lightDistance, spotLight.distance, spotLight.decay );\n\t\t\tdirectLight.visible = true;\n\t\t} else {\n\t\t\tdirectLight.color = vec3( 0.0 );\n\t\t\tdirectLight.visible = false;\n\t\t}\n\t}\n#endif\n#if NUM_RECT_AREA_LIGHTS > 0\n\tstruct RectAreaLight {\n\t\tvec3 color;\n\t\tvec3 position;\n\t\tvec3 halfWidth;\n\t\tvec3 halfHeight;\n\t};\n\tuniform sampler2D ltc_1;\tuniform sampler2D ltc_2;\n\tuniform RectAreaLight rectAreaLights[ NUM_RECT_AREA_LIGHTS ];\n#endif\n#if NUM_HEMI_LIGHTS > 0\n\tstruct HemisphereLight {\n\t\tvec3 direction;\n\t\tvec3 skyColor;\n\t\tvec3 groundColor;\n\t};\n\tuniform HemisphereLight hemisphereLights[ NUM_HEMI_LIGHTS ];\n\tvec3 getHemisphereLightIrradiance( const in HemisphereLight hemiLight, const in GeometricContext geometry ) {\n\t\tfloat dotNL = dot( geometry.normal, hemiLight.direction );\n\t\tfloat hemiDiffuseWeight = 0.5 * dotNL + 0.5;\n\t\tvec3 irradiance = mix( hemiLight.groundColor, hemiLight.skyColor, hemiDiffuseWeight );\n\t\t#ifndef PHYSICALLY_CORRECT_LIGHTS\n\t\t\tirradiance *= PI;\n\t\t#endif\n\t\treturn irradiance;\n\t}\n#endif",
	envmap_physical_pars_fragment = "#if defined( USE_ENVMAP )\n\t#ifdef ENVMAP_MODE_REFRACTION\n\t\tuniform float refractionRatio;\n\t#endif\n\tvec3 getLightProbeIndirectIrradiance( const in GeometricContext geometry, const in int maxMIPLevel ) {\n\t\tvec3 worldNormal = inverseTransformDirection( geometry.normal, viewMatrix );\n\t\t#ifdef ENVMAP_TYPE_CUBE\n\t\t\tvec3 queryVec = vec3( flipEnvMap * worldNormal.x, worldNormal.yz );\n\t\t\t#ifdef TEXTURE_LOD_EXT\n\t\t\t\tvec4 envMapColor = textureCubeLodEXT( envMap, queryVec, float( maxMIPLevel ) );\n\t\t\t#else\n\t\t\t\tvec4 envMapColor = textureCube( envMap, queryVec, float( maxMIPLevel ) );\n\t\t\t#endif\n\t\t\tenvMapColor.rgb = envMapTexelToLinear( envMapColor ).rgb;\n\t\t#elif defined( ENVMAP_TYPE_CUBE_UV )\n\t\t\tvec4 envMapColor = textureCubeUV( envMap, worldNormal, 1.0 );\n\t\t#else\n\t\t\tvec4 envMapColor = vec4( 0.0 );\n\t\t#endif\n\t\treturn PI * envMapColor.rgb * envMapIntensity;\n\t}\n\tfloat getSpecularMIPLevel( const in float roughness, const in int maxMIPLevel ) {\n\t\tfloat maxMIPLevelScalar = float( maxMIPLevel );\n\t\tfloat sigma = PI * roughness * roughness / ( 1.0 + roughness );\n\t\tfloat desiredMIPLevel = maxMIPLevelScalar + log2( sigma );\n\t\treturn clamp( desiredMIPLevel, 0.0, maxMIPLevelScalar );\n\t}\n\tvec3 getLightProbeIndirectRadiance( const in vec3 viewDir, const in vec3 normal, const in float roughness, const in int maxMIPLevel ) {\n\t\t#ifdef ENVMAP_MODE_REFLECTION\n\t\t\tvec3 reflectVec = reflect( -viewDir, normal );\n\t\t\treflectVec = normalize( mix( reflectVec, normal, roughness * roughness) );\n\t\t#else\n\t\t\tvec3 reflectVec = refract( -viewDir, normal, refractionRatio );\n\t\t#endif\n\t\treflectVec = inverseTransformDirection( reflectVec, viewMatrix );\n\t\tfloat specularMIPLevel = getSpecularMIPLevel( roughness, maxMIPLevel );\n\t\t#ifdef ENVMAP_TYPE_CUBE\n\t\t\tvec3 queryReflectVec = vec3( flipEnvMap * reflectVec.x, reflectVec.yz );\n\t\t\t#ifdef TEXTURE_LOD_EXT\n\t\t\t\tvec4 envMapColor = textureCubeLodEXT( envMap, queryReflectVec, specularMIPLevel );\n\t\t\t#else\n\t\t\t\tvec4 envMapColor = textureCube( envMap, queryReflectVec, specularMIPLevel );\n\t\t\t#endif\n\t\t\tenvMapColor.rgb = envMapTexelToLinear( envMapColor ).rgb;\n\t\t#elif defined( ENVMAP_TYPE_CUBE_UV )\n\t\t\tvec4 envMapColor = textureCubeUV( envMap, reflectVec, roughness );\n\t\t#endif\n\t\treturn envMapColor.rgb * envMapIntensity;\n\t}\n#endif",
	lights_toon_fragment = "ToonMaterial material;\nmaterial.diffuseColor = diffuseColor.rgb;",
	lights_toon_pars_fragment = "varying vec3 vViewPosition;\n#ifndef FLAT_SHADED\n\tvarying vec3 vNormal;\n#endif\nstruct ToonMaterial {\n\tvec3 diffuseColor;\n};\nvoid RE_Direct_Toon( const in IncidentLight directLight, const in GeometricContext geometry, const in ToonMaterial material, inout ReflectedLight reflectedLight ) {\n\tvec3 irradiance = getGradientIrradiance( geometry.normal, directLight.direction ) * directLight.color;\n\t#ifndef PHYSICALLY_CORRECT_LIGHTS\n\t\tirradiance *= PI;\n\t#endif\n\treflectedLight.directDiffuse += irradiance * BRDF_Diffuse_Lambert( material.diffuseColor );\n}\nvoid RE_IndirectDiffuse_Toon( const in vec3 irradiance, const in GeometricContext geometry, const in ToonMaterial material, inout ReflectedLight reflectedLight ) {\n\treflectedLight.indirectDiffuse += irradiance * BRDF_Diffuse_Lambert( material.diffuseColor );\n}\n#define RE_Direct\t\t\t\tRE_Direct_Toon\n#define RE_IndirectDiffuse\t\tRE_IndirectDiffuse_Toon\n#define Material_LightProbeLOD( material )\t(0)",
	lights_phong_fragment = "BlinnPhongMaterial material;\nmaterial.diffuseColor = diffuseColor.rgb;\nmaterial.specularColor = specular;\nmaterial.specularShininess = shininess;\nmaterial.specularStrength = specularStrength;",
	lights_phong_pars_fragment = "varying vec3 vViewPosition;\n#ifndef FLAT_SHADED\n\tvarying vec3 vNormal;\n#endif\nstruct BlinnPhongMaterial {\n\tvec3 diffuseColor;\n\tvec3 specularColor;\n\tfloat specularShininess;\n\tfloat specularStrength;\n};\nvoid RE_Direct_BlinnPhong( const in IncidentLight directLight, const in GeometricContext geometry, const in BlinnPhongMaterial material, inout ReflectedLight reflectedLight ) {\n\tfloat dotNL = saturate( dot( geometry.normal, directLight.direction ) );\n\tvec3 irradiance = dotNL * directLight.color;\n\t#ifndef PHYSICALLY_CORRECT_LIGHTS\n\t\tirradiance *= PI;\n\t#endif\n\treflectedLight.directDiffuse += irradiance * BRDF_Diffuse_Lambert( material.diffuseColor );\n\treflectedLight.directSpecular += irradiance * BRDF_Specular_BlinnPhong( directLight, geometry, material.specularColor, material.specularShininess ) * material.specularStrength;\n}\nvoid RE_IndirectDiffuse_BlinnPhong( const in vec3 irradiance, const in GeometricContext geometry, const in BlinnPhongMaterial material, inout ReflectedLight reflectedLight ) {\n\treflectedLight.indirectDiffuse += irradiance * BRDF_Diffuse_Lambert( material.diffuseColor );\n}\n#define RE_Direct\t\t\t\tRE_Direct_BlinnPhong\n#define RE_IndirectDiffuse\t\tRE_IndirectDiffuse_BlinnPhong\n#define Material_LightProbeLOD( material )\t(0)",
	lights_physical_fragment = "PhysicalMaterial material;\nmaterial.diffuseColor = diffuseColor.rgb * ( 1.0 - metalnessFactor );\nvec3 dxy = max( abs( dFdx( geometryNormal ) ), abs( dFdy( geometryNormal ) ) );\nfloat geometryRoughness = max( max( dxy.x, dxy.y ), dxy.z );\nmaterial.specularRoughness = max( roughnessFactor, 0.0525 );material.specularRoughness += geometryRoughness;\nmaterial.specularRoughness = min( material.specularRoughness, 1.0 );\n#ifdef REFLECTIVITY\n\tmaterial.specularColor = mix( vec3( MAXIMUM_SPECULAR_COEFFICIENT * pow2( reflectivity ) ), diffuseColor.rgb, metalnessFactor );\n#else\n\tmaterial.specularColor = mix( vec3( DEFAULT_SPECULAR_COEFFICIENT ), diffuseColor.rgb, metalnessFactor );\n#endif\n#ifdef CLEARCOAT\n\tmaterial.clearcoat = clearcoat;\n\tmaterial.clearcoatRoughness = clearcoatRoughness;\n\t#ifdef USE_CLEARCOATMAP\n\t\tmaterial.clearcoat *= texture2D( clearcoatMap, vUv ).x;\n\t#endif\n\t#ifdef USE_CLEARCOAT_ROUGHNESSMAP\n\t\tmaterial.clearcoatRoughness *= texture2D( clearcoatRoughnessMap, vUv ).y;\n\t#endif\n\tmaterial.clearcoat = saturate( material.clearcoat );\tmaterial.clearcoatRoughness = max( material.clearcoatRoughness, 0.0525 );\n\tmaterial.clearcoatRoughness += geometryRoughness;\n\tmaterial.clearcoatRoughness = min( material.clearcoatRoughness, 1.0 );\n#endif\n#ifdef USE_SHEEN\n\tmaterial.sheenColor = sheen;\n#endif",
	lights_physical_pars_fragment = "struct PhysicalMaterial {\n\tvec3 diffuseColor;\n\tfloat specularRoughness;\n\tvec3 specularColor;\n#ifdef CLEARCOAT\n\tfloat clearcoat;\n\tfloat clearcoatRoughness;\n#endif\n#ifdef USE_SHEEN\n\tvec3 sheenColor;\n#endif\n};\n#define MAXIMUM_SPECULAR_COEFFICIENT 0.16\n#define DEFAULT_SPECULAR_COEFFICIENT 0.04\nfloat clearcoatDHRApprox( const in float roughness, const in float dotNL ) {\n\treturn DEFAULT_SPECULAR_COEFFICIENT + ( 1.0 - DEFAULT_SPECULAR_COEFFICIENT ) * ( pow( 1.0 - dotNL, 5.0 ) * pow( 1.0 - roughness, 2.0 ) );\n}\n#if NUM_RECT_AREA_LIGHTS > 0\n\tvoid RE_Direct_RectArea_Physical( const in RectAreaLight rectAreaLight, const in GeometricContext geometry, const in PhysicalMaterial material, inout ReflectedLight reflectedLight ) {\n\t\tvec3 normal = geometry.normal;\n\t\tvec3 viewDir = geometry.viewDir;\n\t\tvec3 position = geometry.position;\n\t\tvec3 lightPos = rectAreaLight.position;\n\t\tvec3 halfWidth = rectAreaLight.halfWidth;\n\t\tvec3 halfHeight = rectAreaLight.halfHeight;\n\t\tvec3 lightColor = rectAreaLight.color;\n\t\tfloat roughness = material.specularRoughness;\n\t\tvec3 rectCoords[ 4 ];\n\t\trectCoords[ 0 ] = lightPos + halfWidth - halfHeight;\t\trectCoords[ 1 ] = lightPos - halfWidth - halfHeight;\n\t\trectCoords[ 2 ] = lightPos - halfWidth + halfHeight;\n\t\trectCoords[ 3 ] = lightPos + halfWidth + halfHeight;\n\t\tvec2 uv = LTC_Uv( normal, viewDir, roughness );\n\t\tvec4 t1 = texture2D( ltc_1, uv );\n\t\tvec4 t2 = texture2D( ltc_2, uv );\n\t\tmat3 mInv = mat3(\n\t\t\tvec3( t1.x, 0, t1.y ),\n\t\t\tvec3(    0, 1,    0 ),\n\t\t\tvec3( t1.z, 0, t1.w )\n\t\t);\n\t\tvec3 fresnel = ( material.specularColor * t2.x + ( vec3( 1.0 ) - material.specularColor ) * t2.y );\n\t\treflectedLight.directSpecular += lightColor * fresnel * LTC_Evaluate( normal, viewDir, position, mInv, rectCoords );\n\t\treflectedLight.directDiffuse += lightColor * material.diffuseColor * LTC_Evaluate( normal, viewDir, position, mat3( 1.0 ), rectCoords );\n\t}\n#endif\nvoid RE_Direct_Physical( const in IncidentLight directLight, const in GeometricContext geometry, const in PhysicalMaterial material, inout ReflectedLight reflectedLight ) {\n\tfloat dotNL = saturate( dot( geometry.normal, directLight.direction ) );\n\tvec3 irradiance = dotNL * directLight.color;\n\t#ifndef PHYSICALLY_CORRECT_LIGHTS\n\t\tirradiance *= PI;\n\t#endif\n\t#ifdef CLEARCOAT\n\t\tfloat ccDotNL = saturate( dot( geometry.clearcoatNormal, directLight.direction ) );\n\t\tvec3 ccIrradiance = ccDotNL * directLight.color;\n\t\t#ifndef PHYSICALLY_CORRECT_LIGHTS\n\t\t\tccIrradiance *= PI;\n\t\t#endif\n\t\tfloat clearcoatDHR = material.clearcoat * clearcoatDHRApprox( material.clearcoatRoughness, ccDotNL );\n\t\treflectedLight.directSpecular += ccIrradiance * material.clearcoat * BRDF_Specular_GGX( directLight, geometry.viewDir, geometry.clearcoatNormal, vec3( DEFAULT_SPECULAR_COEFFICIENT ), material.clearcoatRoughness );\n\t#else\n\t\tfloat clearcoatDHR = 0.0;\n\t#endif\n\t#ifdef USE_SHEEN\n\t\treflectedLight.directSpecular += ( 1.0 - clearcoatDHR ) * irradiance * BRDF_Specular_Sheen(\n\t\t\tmaterial.specularRoughness,\n\t\t\tdirectLight.direction,\n\t\t\tgeometry,\n\t\t\tmaterial.sheenColor\n\t\t);\n\t#else\n\t\treflectedLight.directSpecular += ( 1.0 - clearcoatDHR ) * irradiance * BRDF_Specular_GGX( directLight, geometry.viewDir, geometry.normal, material.specularColor, material.specularRoughness);\n\t#endif\n\treflectedLight.directDiffuse += ( 1.0 - clearcoatDHR ) * irradiance * BRDF_Diffuse_Lambert( material.diffuseColor );\n}\nvoid RE_IndirectDiffuse_Physical( const in vec3 irradiance, const in GeometricContext geometry, const in PhysicalMaterial material, inout ReflectedLight reflectedLight ) {\n\treflectedLight.indirectDiffuse += irradiance * BRDF_Diffuse_Lambert( material.diffuseColor );\n}\nvoid RE_IndirectSpecular_Physical( const in vec3 radiance, const in vec3 irradiance, const in vec3 clearcoatRadiance, const in GeometricContext geometry, const in PhysicalMaterial material, inout ReflectedLight reflectedLight) {\n\t#ifdef CLEARCOAT\n\t\tfloat ccDotNV = saturate( dot( geometry.clearcoatNormal, geometry.viewDir ) );\n\t\treflectedLight.indirectSpecular += clearcoatRadiance * material.clearcoat * BRDF_Specular_GGX_Environment( geometry.viewDir, geometry.clearcoatNormal, vec3( DEFAULT_SPECULAR_COEFFICIENT ), material.clearcoatRoughness );\n\t\tfloat ccDotNL = ccDotNV;\n\t\tfloat clearcoatDHR = material.clearcoat * clearcoatDHRApprox( material.clearcoatRoughness, ccDotNL );\n\t#else\n\t\tfloat clearcoatDHR = 0.0;\n\t#endif\n\tfloat clearcoatInv = 1.0 - clearcoatDHR;\n\tvec3 singleScattering = vec3( 0.0 );\n\tvec3 multiScattering = vec3( 0.0 );\n\tvec3 cosineWeightedIrradiance = irradiance * RECIPROCAL_PI;\n\tBRDF_Specular_Multiscattering_Environment( geometry, material.specularColor, material.specularRoughness, singleScattering, multiScattering );\n\tvec3 diffuse = material.diffuseColor * ( 1.0 - ( singleScattering + multiScattering ) );\n\treflectedLight.indirectSpecular += clearcoatInv * radiance * singleScattering;\n\treflectedLight.indirectSpecular += multiScattering * cosineWeightedIrradiance;\n\treflectedLight.indirectDiffuse += diffuse * cosineWeightedIrradiance;\n}\n#define RE_Direct\t\t\t\tRE_Direct_Physical\n#define RE_Direct_RectArea\t\tRE_Direct_RectArea_Physical\n#define RE_IndirectDiffuse\t\tRE_IndirectDiffuse_Physical\n#define RE_IndirectSpecular\t\tRE_IndirectSpecular_Physical\nfloat computeSpecularOcclusion( const in float dotNV, const in float ambientOcclusion, const in float roughness ) {\n\treturn saturate( pow( dotNV + ambientOcclusion, exp2( - 16.0 * roughness - 1.0 ) ) - 1.0 + ambientOcclusion );\n}",
	lights_fragment_begin = "\nGeometricContext geometry;\ngeometry.position = - vViewPosition;\ngeometry.normal = normal;\ngeometry.viewDir = ( isOrthographic ) ? vec3( 0, 0, 1 ) : normalize( vViewPosition );\n#ifdef CLEARCOAT\n\tgeometry.clearcoatNormal = clearcoatNormal;\n#endif\nIncidentLight directLight;\n#if ( NUM_POINT_LIGHTS > 0 ) && defined( RE_Direct )\n\tPointLight pointLight;\n\t#if defined( USE_SHADOWMAP ) && NUM_POINT_LIGHT_SHADOWS > 0\n\tPointLightShadow pointLightShadow;\n\t#endif\n\t#pragma unroll_loop_start\n\tfor ( int i = 0; i < NUM_POINT_LIGHTS; i ++ ) {\n\t\tpointLight = pointLights[ i ];\n\t\tgetPointDirectLightIrradiance( pointLight, geometry, directLight );\n\t\t#if defined( USE_SHADOWMAP ) && ( UNROLLED_LOOP_INDEX < NUM_POINT_LIGHT_SHADOWS )\n\t\tpointLightShadow = pointLightShadows[ i ];\n\t\tdirectLight.color *= all( bvec2( directLight.visible, receiveShadow ) ) ? getPointShadow( pointShadowMap[ i ], pointLightShadow.shadowMapSize, pointLightShadow.shadowBias, pointLightShadow.shadowRadius, vPointShadowCoord[ i ], pointLightShadow.shadowCameraNear, pointLightShadow.shadowCameraFar ) : 1.0;\n\t\t#endif\n\t\tRE_Direct( directLight, geometry, material, reflectedLight );\n\t}\n\t#pragma unroll_loop_end\n#endif\n#if ( NUM_SPOT_LIGHTS > 0 ) && defined( RE_Direct )\n\tSpotLight spotLight;\n\t#if defined( USE_SHADOWMAP ) && NUM_SPOT_LIGHT_SHADOWS > 0\n\tSpotLightShadow spotLightShadow;\n\t#endif\n\t#pragma unroll_loop_start\n\tfor ( int i = 0; i < NUM_SPOT_LIGHTS; i ++ ) {\n\t\tspotLight = spotLights[ i ];\n\t\tgetSpotDirectLightIrradiance( spotLight, geometry, directLight );\n\t\t#if defined( USE_SHADOWMAP ) && ( UNROLLED_LOOP_INDEX < NUM_SPOT_LIGHT_SHADOWS )\n\t\tspotLightShadow = spotLightShadows[ i ];\n\t\tdirectLight.color *= all( bvec2( directLight.visible, receiveShadow ) ) ? getShadow( spotShadowMap[ i ], spotLightShadow.shadowMapSize, spotLightShadow.shadowBias, spotLightShadow.shadowRadius, vSpotShadowCoord[ i ] ) : 1.0;\n\t\t#endif\n\t\tRE_Direct( directLight, geometry, material, reflectedLight );\n\t}\n\t#pragma unroll_loop_end\n#endif\n#if ( NUM_DIR_LIGHTS > 0 ) && defined( RE_Direct )\n\tDirectionalLight directionalLight;\n\t#if defined( USE_SHADOWMAP ) && NUM_DIR_LIGHT_SHADOWS > 0\n\tDirectionalLightShadow directionalLightShadow;\n\t#endif\n\t#pragma unroll_loop_start\n\tfor ( int i = 0; i < NUM_DIR_LIGHTS; i ++ ) {\n\t\tdirectionalLight = directionalLights[ i ];\n\t\tgetDirectionalDirectLightIrradiance( directionalLight, geometry, directLight );\n\t\t#if defined( USE_SHADOWMAP ) && ( UNROLLED_LOOP_INDEX < NUM_DIR_LIGHT_SHADOWS )\n\t\tdirectionalLightShadow = directionalLightShadows[ i ];\n\t\tdirectLight.color *= all( bvec2( directLight.visible, receiveShadow ) ) ? getShadow( directionalShadowMap[ i ], directionalLightShadow.shadowMapSize, directionalLightShadow.shadowBias, directionalLightShadow.shadowRadius, vDirectionalShadowCoord[ i ] ) : 1.0;\n\t\t#endif\n\t\tRE_Direct( directLight, geometry, material, reflectedLight );\n\t}\n\t#pragma unroll_loop_end\n#endif\n#if ( NUM_RECT_AREA_LIGHTS > 0 ) && defined( RE_Direct_RectArea )\n\tRectAreaLight rectAreaLight;\n\t#pragma unroll_loop_start\n\tfor ( int i = 0; i < NUM_RECT_AREA_LIGHTS; i ++ ) {\n\t\trectAreaLight = rectAreaLights[ i ];\n\t\tRE_Direct_RectArea( rectAreaLight, geometry, material, reflectedLight );\n\t}\n\t#pragma unroll_loop_end\n#endif\n#if defined( RE_IndirectDiffuse )\n\tvec3 iblIrradiance = vec3( 0.0 );\n\tvec3 irradiance = getAmbientLightIrradiance( ambientLightColor );\n\tirradiance += getLightProbeIrradiance( lightProbe, geometry );\n\t#if ( NUM_HEMI_LIGHTS > 0 )\n\t\t#pragma unroll_loop_start\n\t\tfor ( int i = 0; i < NUM_HEMI_LIGHTS; i ++ ) {\n\t\t\tirradiance += getHemisphereLightIrradiance( hemisphereLights[ i ], geometry );\n\t\t}\n\t\t#pragma unroll_loop_end\n\t#endif\n#endif\n#if defined( RE_IndirectSpecular )\n\tvec3 radiance = vec3( 0.0 );\n\tvec3 clearcoatRadiance = vec3( 0.0 );\n#endif",
	lights_fragment_maps = "#if defined( RE_IndirectDiffuse )\n\t#ifdef USE_LIGHTMAP\n\t\tvec4 lightMapTexel= texture2D( lightMap, vUv2 );\n\t\tvec3 lightMapIrradiance = lightMapTexelToLinear( lightMapTexel ).rgb * lightMapIntensity;\n\t\t#ifndef PHYSICALLY_CORRECT_LIGHTS\n\t\t\tlightMapIrradiance *= PI;\n\t\t#endif\n\t\tirradiance += lightMapIrradiance;\n\t#endif\n\t#if defined( USE_ENVMAP ) && defined( STANDARD ) && defined( ENVMAP_TYPE_CUBE_UV )\n\t\tiblIrradiance += getLightProbeIndirectIrradiance( geometry, maxMipLevel );\n\t#endif\n#endif\n#if defined( USE_ENVMAP ) && defined( RE_IndirectSpecular )\n\tradiance += getLightProbeIndirectRadiance( geometry.viewDir, geometry.normal, material.specularRoughness, maxMipLevel );\n\t#ifdef CLEARCOAT\n\t\tclearcoatRadiance += getLightProbeIndirectRadiance( geometry.viewDir, geometry.clearcoatNormal, material.clearcoatRoughness, maxMipLevel );\n\t#endif\n#endif",
	lights_fragment_end = "#if defined( RE_IndirectDiffuse )\n\tRE_IndirectDiffuse( irradiance, geometry, material, reflectedLight );\n#endif\n#if defined( RE_IndirectSpecular )\n\tRE_IndirectSpecular( radiance, iblIrradiance, clearcoatRadiance, geometry, material, reflectedLight );\n#endif",
	logdepthbuf_fragment = "#if defined( USE_LOGDEPTHBUF ) && defined( USE_LOGDEPTHBUF_EXT )\n\tgl_FragDepthEXT = vIsPerspective == 0.0 ? gl_FragCoord.z : log2( vFragDepth ) * logDepthBufFC * 0.5;\n#endif",
	logdepthbuf_pars_fragment = "#if defined( USE_LOGDEPTHBUF ) && defined( USE_LOGDEPTHBUF_EXT )\n\tuniform float logDepthBufFC;\n\tvarying float vFragDepth;\n\tvarying float vIsPerspective;\n#endif",
	logdepthbuf_pars_vertex = "#ifdef USE_LOGDEPTHBUF\n\t#ifdef USE_LOGDEPTHBUF_EXT\n\t\tvarying float vFragDepth;\n\t\tvarying float vIsPerspective;\n\t#else\n\t\tuniform float logDepthBufFC;\n\t#endif\n#endif",
	logdepthbuf_vertex = "#ifdef USE_LOGDEPTHBUF\n\t#ifdef USE_LOGDEPTHBUF_EXT\n\t\tvFragDepth = 1.0 + gl_Position.w;\n\t\tvIsPerspective = float( isPerspectiveMatrix( projectionMatrix ) );\n\t#else\n\t\tif ( isPerspectiveMatrix( projectionMatrix ) ) {\n\t\t\tgl_Position.z = log2( max( EPSILON, gl_Position.w + 1.0 ) ) * logDepthBufFC - 1.0;\n\t\t\tgl_Position.z *= gl_Position.w;\n\t\t}\n\t#endif\n#endif",
	map_fragment = "#ifdef USE_MAP\n\tvec4 texelColor = texture2D( map, vUv );\n\ttexelColor = mapTexelToLinear( texelColor );\n\tdiffuseColor *= texelColor;\n#endif",
	map_pars_fragment = "#ifdef USE_MAP\n\tuniform sampler2D map;\n#endif",
	map_particle_fragment = "#if defined( USE_MAP ) || defined( USE_ALPHAMAP )\n\tvec2 uv = ( uvTransform * vec3( gl_PointCoord.x, 1.0 - gl_PointCoord.y, 1 ) ).xy;\n#endif\n#ifdef USE_MAP\n\tvec4 mapTexel = texture2D( map, uv );\n\tdiffuseColor *= mapTexelToLinear( mapTexel );\n#endif\n#ifdef USE_ALPHAMAP\n\tdiffuseColor.a *= texture2D( alphaMap, uv ).g;\n#endif",
	map_particle_pars_fragment = "#if defined( USE_MAP ) || defined( USE_ALPHAMAP )\n\tuniform mat3 uvTransform;\n#endif\n#ifdef USE_MAP\n\tuniform sampler2D map;\n#endif\n#ifdef USE_ALPHAMAP\n\tuniform sampler2D alphaMap;\n#endif",
	metalnessmap_fragment = "float metalnessFactor = metalness;\n#ifdef USE_METALNESSMAP\n\tvec4 texelMetalness = texture2D( metalnessMap, vUv );\n\tmetalnessFactor *= texelMetalness.b;\n#endif",
	metalnessmap_pars_fragment = "#ifdef USE_METALNESSMAP\n\tuniform sampler2D metalnessMap;\n#endif",
	morphnormal_vertex = "#ifdef USE_MORPHNORMALS\n\tobjectNormal *= morphTargetBaseInfluence;\n\tobjectNormal += morphNormal0 * morphTargetInfluences[ 0 ];\n\tobjectNormal += morphNormal1 * morphTargetInfluences[ 1 ];\n\tobjectNormal += morphNormal2 * morphTargetInfluences[ 2 ];\n\tobjectNormal += morphNormal3 * morphTargetInfluences[ 3 ];\n#endif",
	morphtarget_pars_vertex = "#ifdef USE_MORPHTARGETS\n\tuniform float morphTargetBaseInfluence;\n\t#ifndef USE_MORPHNORMALS\n\t\tuniform float morphTargetInfluences[ 8 ];\n\t#else\n\t\tuniform float morphTargetInfluences[ 4 ];\n\t#endif\n#endif",
	morphtarget_vertex = "#ifdef USE_MORPHTARGETS\n\ttransformed *= morphTargetBaseInfluence;\n\ttransformed += morphTarget0 * morphTargetInfluences[ 0 ];\n\ttransformed += morphTarget1 * morphTargetInfluences[ 1 ];\n\ttransformed += morphTarget2 * morphTargetInfluences[ 2 ];\n\ttransformed += morphTarget3 * morphTargetInfluences[ 3 ];\n\t#ifndef USE_MORPHNORMALS\n\t\ttransformed += morphTarget4 * morphTargetInfluences[ 4 ];\n\t\ttransformed += morphTarget5 * morphTargetInfluences[ 5 ];\n\t\ttransformed += morphTarget6 * morphTargetInfluences[ 6 ];\n\t\ttransformed += morphTarget7 * morphTargetInfluences[ 7 ];\n\t#endif\n#endif",
	normal_fragment_begin = "float faceDirection = gl_FrontFacing ? 1.0 : - 1.0;\n#ifdef FLAT_SHADED\n\tvec3 fdx = vec3( dFdx( vViewPosition.x ), dFdx( vViewPosition.y ), dFdx( vViewPosition.z ) );\n\tvec3 fdy = vec3( dFdy( vViewPosition.x ), dFdy( vViewPosition.y ), dFdy( vViewPosition.z ) );\n\tvec3 normal = normalize( cross( fdx, fdy ) );\n#else\n\tvec3 normal = normalize( vNormal );\n\t#ifdef DOUBLE_SIDED\n\t\tnormal = normal * faceDirection;\n\t#endif\n\t#ifdef USE_TANGENT\n\t\tvec3 tangent = normalize( vTangent );\n\t\tvec3 bitangent = normalize( vBitangent );\n\t\t#ifdef DOUBLE_SIDED\n\t\t\ttangent = tangent * faceDirection;\n\t\t\tbitangent = bitangent * faceDirection;\n\t\t#endif\n\t\t#if defined( TANGENTSPACE_NORMALMAP ) || defined( USE_CLEARCOAT_NORMALMAP )\n\t\t\tmat3 vTBN = mat3( tangent, bitangent, normal );\n\t\t#endif\n\t#endif\n#endif\nvec3 geometryNormal = normal;",
	normal_fragment_maps = "#ifdef OBJECTSPACE_NORMALMAP\n\tnormal = texture2D( normalMap, vUv ).xyz * 2.0 - 1.0;\n\t#ifdef FLIP_SIDED\n\t\tnormal = - normal;\n\t#endif\n\t#ifdef DOUBLE_SIDED\n\t\tnormal = normal * faceDirection;\n\t#endif\n\tnormal = normalize( normalMatrix * normal );\n#elif defined( TANGENTSPACE_NORMALMAP )\n\tvec3 mapN = texture2D( normalMap, vUv ).xyz * 2.0 - 1.0;\n\tmapN.xy *= normalScale;\n\t#ifdef USE_TANGENT\n\t\tnormal = normalize( vTBN * mapN );\n\t#else\n\t\tnormal = perturbNormal2Arb( -vViewPosition, normal, mapN, faceDirection );\n\t#endif\n#elif defined( USE_BUMPMAP )\n\tnormal = perturbNormalArb( -vViewPosition, normal, dHdxy_fwd(), faceDirection );\n#endif",
	normalmap_pars_fragment = "#ifdef USE_NORMALMAP\n\tuniform sampler2D normalMap;\n\tuniform vec2 normalScale;\n#endif\n#ifdef OBJECTSPACE_NORMALMAP\n\tuniform mat3 normalMatrix;\n#endif\n#if ! defined ( USE_TANGENT ) && ( defined ( TANGENTSPACE_NORMALMAP ) || defined ( USE_CLEARCOAT_NORMALMAP ) )\n\tvec3 perturbNormal2Arb( vec3 eye_pos, vec3 surf_norm, vec3 mapN, float faceDirection ) {\n\t\tvec3 q0 = vec3( dFdx( eye_pos.x ), dFdx( eye_pos.y ), dFdx( eye_pos.z ) );\n\t\tvec3 q1 = vec3( dFdy( eye_pos.x ), dFdy( eye_pos.y ), dFdy( eye_pos.z ) );\n\t\tvec2 st0 = dFdx( vUv.st );\n\t\tvec2 st1 = dFdy( vUv.st );\n\t\tvec3 N = surf_norm;\n\t\tvec3 q1perp = cross( q1, N );\n\t\tvec3 q0perp = cross( N, q0 );\n\t\tvec3 T = q1perp * st0.x + q0perp * st1.x;\n\t\tvec3 B = q1perp * st0.y + q0perp * st1.y;\n\t\tfloat det = max( dot( T, T ), dot( B, B ) );\n\t\tfloat scale = ( det == 0.0 ) ? 0.0 : faceDirection * inversesqrt( det );\n\t\treturn normalize( T * ( mapN.x * scale ) + B * ( mapN.y * scale ) + N * mapN.z );\n\t}\n#endif",
	clearcoat_normal_fragment_begin = "#ifdef CLEARCOAT\n\tvec3 clearcoatNormal = geometryNormal;\n#endif",
	clearcoat_normal_fragment_maps = "#ifdef USE_CLEARCOAT_NORMALMAP\n\tvec3 clearcoatMapN = texture2D( clearcoatNormalMap, vUv ).xyz * 2.0 - 1.0;\n\tclearcoatMapN.xy *= clearcoatNormalScale;\n\t#ifdef USE_TANGENT\n\t\tclearcoatNormal = normalize( vTBN * clearcoatMapN );\n\t#else\n\t\tclearcoatNormal = perturbNormal2Arb( - vViewPosition, clearcoatNormal, clearcoatMapN, faceDirection );\n\t#endif\n#endif",
	clearcoat_pars_fragment = "#ifdef USE_CLEARCOATMAP\n\tuniform sampler2D clearcoatMap;\n#endif\n#ifdef USE_CLEARCOAT_ROUGHNESSMAP\n\tuniform sampler2D clearcoatRoughnessMap;\n#endif\n#ifdef USE_CLEARCOAT_NORMALMAP\n\tuniform sampler2D clearcoatNormalMap;\n\tuniform vec2 clearcoatNormalScale;\n#endif",
	packing = "vec3 packNormalToRGB( const in vec3 normal ) {\n\treturn normalize( normal ) * 0.5 + 0.5;\n}\nvec3 unpackRGBToNormal( const in vec3 rgb ) {\n\treturn 2.0 * rgb.xyz - 1.0;\n}\nconst float PackUpscale = 256. / 255.;const float UnpackDownscale = 255. / 256.;\nconst vec3 PackFactors = vec3( 256. * 256. * 256., 256. * 256., 256. );\nconst vec4 UnpackFactors = UnpackDownscale / vec4( PackFactors, 1. );\nconst float ShiftRight8 = 1. / 256.;\nvec4 packDepthToRGBA( const in float v ) {\n\tvec4 r = vec4( fract( v * PackFactors ), v );\n\tr.yzw -= r.xyz * ShiftRight8;\treturn r * PackUpscale;\n}\nfloat unpackRGBAToDepth( const in vec4 v ) {\n\treturn dot( v, UnpackFactors );\n}\nvec4 pack2HalfToRGBA( vec2 v ) {\n\tvec4 r = vec4( v.x, fract( v.x * 255.0 ), v.y, fract( v.y * 255.0 ));\n\treturn vec4( r.x - r.y / 255.0, r.y, r.z - r.w / 255.0, r.w);\n}\nvec2 unpackRGBATo2Half( vec4 v ) {\n\treturn vec2( v.x + ( v.y / 255.0 ), v.z + ( v.w / 255.0 ) );\n}\nfloat viewZToOrthographicDepth( const in float viewZ, const in float near, const in float far ) {\n\treturn ( viewZ + near ) / ( near - far );\n}\nfloat orthographicDepthToViewZ( const in float linearClipZ, const in float near, const in float far ) {\n\treturn linearClipZ * ( near - far ) - near;\n}\nfloat viewZToPerspectiveDepth( const in float viewZ, const in float near, const in float far ) {\n\treturn (( near + viewZ ) * far ) / (( far - near ) * viewZ );\n}\nfloat perspectiveDepthToViewZ( const in float invClipZ, const in float near, const in float far ) {\n\treturn ( near * far ) / ( ( far - near ) * invClipZ - far );\n}",
	premultiplied_alpha_fragment = "#ifdef PREMULTIPLIED_ALPHA\n\tgl_FragColor.rgb *= gl_FragColor.a;\n#endif",
	project_vertex = "vec4 mvPosition = vec4( transformed, 1.0 );\n#ifdef USE_INSTANCING\n\tmvPosition = instanceMatrix * mvPosition;\n#endif\nmvPosition = modelViewMatrix * mvPosition;\ngl_Position = projectionMatrix * mvPosition;",
	dithering_fragment = "#ifdef DITHERING\n\tgl_FragColor.rgb = dithering( gl_FragColor.rgb );\n#endif",
	dithering_pars_fragment = "#ifdef DITHERING\n\tvec3 dithering( vec3 color ) {\n\t\tfloat grid_position = rand( gl_FragCoord.xy );\n\t\tvec3 dither_shift_RGB = vec3( 0.25 / 255.0, -0.25 / 255.0, 0.25 / 255.0 );\n\t\tdither_shift_RGB = mix( 2.0 * dither_shift_RGB, -2.0 * dither_shift_RGB, grid_position );\n\t\treturn color + dither_shift_RGB;\n\t}\n#endif",
	roughnessmap_fragment = "float roughnessFactor = roughness;\n#ifdef USE_ROUGHNESSMAP\n\tvec4 texelRoughness = texture2D( roughnessMap, vUv );\n\troughnessFactor *= texelRoughness.g;\n#endif",
	roughnessmap_pars_fragment = "#ifdef USE_ROUGHNESSMAP\n\tuniform sampler2D roughnessMap;\n#endif",
	shadowmap_pars_fragment = "#ifdef USE_SHADOWMAP\n\t#if NUM_DIR_LIGHT_SHADOWS > 0\n\t\tuniform sampler2D directionalShadowMap[ NUM_DIR_LIGHT_SHADOWS ];\n\t\tvarying vec4 vDirectionalShadowCoord[ NUM_DIR_LIGHT_SHADOWS ];\n\t\tstruct DirectionalLightShadow {\n\t\t\tfloat shadowBias;\n\t\t\tfloat shadowNormalBias;\n\t\t\tfloat shadowRadius;\n\t\t\tvec2 shadowMapSize;\n\t\t};\n\t\tuniform DirectionalLightShadow directionalLightShadows[ NUM_DIR_LIGHT_SHADOWS ];\n\t#endif\n\t#if NUM_SPOT_LIGHT_SHADOWS > 0\n\t\tuniform sampler2D spotShadowMap[ NUM_SPOT_LIGHT_SHADOWS ];\n\t\tvarying vec4 vSpotShadowCoord[ NUM_SPOT_LIGHT_SHADOWS ];\n\t\tstruct SpotLightShadow {\n\t\t\tfloat shadowBias;\n\t\t\tfloat shadowNormalBias;\n\t\t\tfloat shadowRadius;\n\t\t\tvec2 shadowMapSize;\n\t\t};\n\t\tuniform SpotLightShadow spotLightShadows[ NUM_SPOT_LIGHT_SHADOWS ];\n\t#endif\n\t#if NUM_POINT_LIGHT_SHADOWS > 0\n\t\tuniform sampler2D pointShadowMap[ NUM_POINT_LIGHT_SHADOWS ];\n\t\tvarying vec4 vPointShadowCoord[ NUM_POINT_LIGHT_SHADOWS ];\n\t\tstruct PointLightShadow {\n\t\t\tfloat shadowBias;\n\t\t\tfloat shadowNormalBias;\n\t\t\tfloat shadowRadius;\n\t\t\tvec2 shadowMapSize;\n\t\t\tfloat shadowCameraNear;\n\t\t\tfloat shadowCameraFar;\n\t\t};\n\t\tuniform PointLightShadow pointLightShadows[ NUM_POINT_LIGHT_SHADOWS ];\n\t#endif\n\tfloat texture2DCompare( sampler2D depths, vec2 uv, float compare ) {\n\t\treturn step( compare, unpackRGBAToDepth( texture2D( depths, uv ) ) );\n\t}\n\tvec2 texture2DDistribution( sampler2D shadow, vec2 uv ) {\n\t\treturn unpackRGBATo2Half( texture2D( shadow, uv ) );\n\t}\n\tfloat VSMShadow (sampler2D shadow, vec2 uv, float compare ){\n\t\tfloat occlusion = 1.0;\n\t\tvec2 distribution = texture2DDistribution( shadow, uv );\n\t\tfloat hard_shadow = step( compare , distribution.x );\n\t\tif (hard_shadow != 1.0 ) {\n\t\t\tfloat distance = compare - distribution.x ;\n\t\t\tfloat variance = max( 0.00000, distribution.y * distribution.y );\n\t\t\tfloat softness_probability = variance / (variance + distance * distance );\t\t\tsoftness_probability = clamp( ( softness_probability - 0.3 ) / ( 0.95 - 0.3 ), 0.0, 1.0 );\t\t\tocclusion = clamp( max( hard_shadow, softness_probability ), 0.0, 1.0 );\n\t\t}\n\t\treturn occlusion;\n\t}\n\tfloat getShadow( sampler2D shadowMap, vec2 shadowMapSize, float shadowBias, float shadowRadius, vec4 shadowCoord ) {\n\t\tfloat shadow = 1.0;\n\t\tshadowCoord.xyz /= shadowCoord.w;\n\t\tshadowCoord.z += shadowBias;\n\t\tbvec4 inFrustumVec = bvec4 ( shadowCoord.x >= 0.0, shadowCoord.x <= 1.0, shadowCoord.y >= 0.0, shadowCoord.y <= 1.0 );\n\t\tbool inFrustum = all( inFrustumVec );\n\t\tbvec2 frustumTestVec = bvec2( inFrustum, shadowCoord.z <= 1.0 );\n\t\tbool frustumTest = all( frustumTestVec );\n\t\tif ( frustumTest ) {\n\t\t#if defined( SHADOWMAP_TYPE_PCF )\n\t\t\tvec2 texelSize = vec2( 1.0 ) / shadowMapSize;\n\t\t\tfloat dx0 = - texelSize.x * shadowRadius;\n\t\t\tfloat dy0 = - texelSize.y * shadowRadius;\n\t\t\tfloat dx1 = + texelSize.x * shadowRadius;\n\t\t\tfloat dy1 = + texelSize.y * shadowRadius;\n\t\t\tfloat dx2 = dx0 / 2.0;\n\t\t\tfloat dy2 = dy0 / 2.0;\n\t\t\tfloat dx3 = dx1 / 2.0;\n\t\t\tfloat dy3 = dy1 / 2.0;\n\t\t\tshadow = (\n\t\t\t\ttexture2DCompare( shadowMap, shadowCoord.xy + vec2( dx0, dy0 ), shadowCoord.z ) +\n\t\t\t\ttexture2DCompare( shadowMap, shadowCoord.xy + vec2( 0.0, dy0 ), shadowCoord.z ) +\n\t\t\t\ttexture2DCompare( shadowMap, shadowCoord.xy + vec2( dx1, dy0 ), shadowCoord.z ) +\n\t\t\t\ttexture2DCompare( shadowMap, shadowCoord.xy + vec2( dx2, dy2 ), shadowCoord.z ) +\n\t\t\t\ttexture2DCompare( shadowMap, shadowCoord.xy + vec2( 0.0, dy2 ), shadowCoord.z ) +\n\t\t\t\ttexture2DCompare( shadowMap, shadowCoord.xy + vec2( dx3, dy2 ), shadowCoord.z ) +\n\t\t\t\ttexture2DCompare( shadowMap, shadowCoord.xy + vec2( dx0, 0.0 ), shadowCoord.z ) +\n\t\t\t\ttexture2DCompare( shadowMap, shadowCoord.xy + vec2( dx2, 0.0 ), shadowCoord.z ) +\n\t\t\t\ttexture2DCompare( shadowMap, shadowCoord.xy, shadowCoord.z ) +\n\t\t\t\ttexture2DCompare( shadowMap, shadowCoord.xy + vec2( dx3, 0.0 ), shadowCoord.z ) +\n\t\t\t\ttexture2DCompare( shadowMap, shadowCoord.xy + vec2( dx1, 0.0 ), shadowCoord.z ) +\n\t\t\t\ttexture2DCompare( shadowMap, shadowCoord.xy + vec2( dx2, dy3 ), shadowCoord.z ) +\n\t\t\t\ttexture2DCompare( shadowMap, shadowCoord.xy + vec2( 0.0, dy3 ), shadowCoord.z ) +\n\t\t\t\ttexture2DCompare( shadowMap, shadowCoord.xy + vec2( dx3, dy3 ), shadowCoord.z ) +\n\t\t\t\ttexture2DCompare( shadowMap, shadowCoord.xy + vec2( dx0, dy1 ), shadowCoord.z ) +\n\t\t\t\ttexture2DCompare( shadowMap, shadowCoord.xy + vec2( 0.0, dy1 ), shadowCoord.z ) +\n\t\t\t\ttexture2DCompare( shadowMap, shadowCoord.xy + vec2( dx1, dy1 ), shadowCoord.z )\n\t\t\t) * ( 1.0 / 17.0 );\n\t\t#elif defined( SHADOWMAP_TYPE_PCF_SOFT )\n\t\t\tvec2 texelSize = vec2( 1.0 ) / shadowMapSize;\n\t\t\tfloat dx = texelSize.x;\n\t\t\tfloat dy = texelSize.y;\n\t\t\tvec2 uv = shadowCoord.xy;\n\t\t\tvec2 f = fract( uv * shadowMapSize + 0.5 );\n\t\t\tuv -= f * texelSize;\n\t\t\tshadow = (\n\t\t\t\ttexture2DCompare( shadowMap, uv, shadowCoord.z ) +\n\t\t\t\ttexture2DCompare( shadowMap, uv + vec2( dx, 0.0 ), shadowCoord.z ) +\n\t\t\t\ttexture2DCompare( shadowMap, uv + vec2( 0.0, dy ), shadowCoord.z ) +\n\t\t\t\ttexture2DCompare( shadowMap, uv + texelSize, shadowCoord.z ) +\n\t\t\t\tmix( texture2DCompare( shadowMap, uv + vec2( -dx, 0.0 ), shadowCoord.z ), \n\t\t\t\t\t texture2DCompare( shadowMap, uv + vec2( 2.0 * dx, 0.0 ), shadowCoord.z ),\n\t\t\t\t\t f.x ) +\n\t\t\t\tmix( texture2DCompare( shadowMap, uv + vec2( -dx, dy ), shadowCoord.z ), \n\t\t\t\t\t texture2DCompare( shadowMap, uv + vec2( 2.0 * dx, dy ), shadowCoord.z ),\n\t\t\t\t\t f.x ) +\n\t\t\t\tmix( texture2DCompare( shadowMap, uv + vec2( 0.0, -dy ), shadowCoord.z ), \n\t\t\t\t\t texture2DCompare( shadowMap, uv + vec2( 0.0, 2.0 * dy ), shadowCoord.z ),\n\t\t\t\t\t f.y ) +\n\t\t\t\tmix( texture2DCompare( shadowMap, uv + vec2( dx, -dy ), shadowCoord.z ), \n\t\t\t\t\t texture2DCompare( shadowMap, uv + vec2( dx, 2.0 * dy ), shadowCoord.z ),\n\t\t\t\t\t f.y ) +\n\t\t\t\tmix( mix( texture2DCompare( shadowMap, uv + vec2( -dx, -dy ), shadowCoord.z ), \n\t\t\t\t\t\t  texture2DCompare( shadowMap, uv + vec2( 2.0 * dx, -dy ), shadowCoord.z ),\n\t\t\t\t\t\t  f.x ),\n\t\t\t\t\t mix( texture2DCompare( shadowMap, uv + vec2( -dx, 2.0 * dy ), shadowCoord.z ), \n\t\t\t\t\t\t  texture2DCompare( shadowMap, uv + vec2( 2.0 * dx, 2.0 * dy ), shadowCoord.z ),\n\t\t\t\t\t\t  f.x ),\n\t\t\t\t\t f.y )\n\t\t\t) * ( 1.0 / 9.0 );\n\t\t#elif defined( SHADOWMAP_TYPE_VSM )\n\t\t\tshadow = VSMShadow( shadowMap, shadowCoord.xy, shadowCoord.z );\n\t\t#else\n\t\t\tshadow = texture2DCompare( shadowMap, shadowCoord.xy, shadowCoord.z );\n\t\t#endif\n\t\t}\n\t\treturn shadow;\n\t}\n\tvec2 cubeToUV( vec3 v, float texelSizeY ) {\n\t\tvec3 absV = abs( v );\n\t\tfloat scaleToCube = 1.0 / max( absV.x, max( absV.y, absV.z ) );\n\t\tabsV *= scaleToCube;\n\t\tv *= scaleToCube * ( 1.0 - 2.0 * texelSizeY );\n\t\tvec2 planar = v.xy;\n\t\tfloat almostATexel = 1.5 * texelSizeY;\n\t\tfloat almostOne = 1.0 - almostATexel;\n\t\tif ( absV.z >= almostOne ) {\n\t\t\tif ( v.z > 0.0 )\n\t\t\t\tplanar.x = 4.0 - v.x;\n\t\t} else if ( absV.x >= almostOne ) {\n\t\t\tfloat signX = sign( v.x );\n\t\t\tplanar.x = v.z * signX + 2.0 * signX;\n\t\t} else if ( absV.y >= almostOne ) {\n\t\t\tfloat signY = sign( v.y );\n\t\t\tplanar.x = v.x + 2.0 * signY + 2.0;\n\t\t\tplanar.y = v.z * signY - 2.0;\n\t\t}\n\t\treturn vec2( 0.125, 0.25 ) * planar + vec2( 0.375, 0.75 );\n\t}\n\tfloat getPointShadow( sampler2D shadowMap, vec2 shadowMapSize, float shadowBias, float shadowRadius, vec4 shadowCoord, float shadowCameraNear, float shadowCameraFar ) {\n\t\tvec2 texelSize = vec2( 1.0 ) / ( shadowMapSize * vec2( 4.0, 2.0 ) );\n\t\tvec3 lightToPosition = shadowCoord.xyz;\n\t\tfloat dp = ( length( lightToPosition ) - shadowCameraNear ) / ( shadowCameraFar - shadowCameraNear );\t\tdp += shadowBias;\n\t\tvec3 bd3D = normalize( lightToPosition );\n\t\t#if defined( SHADOWMAP_TYPE_PCF ) || defined( SHADOWMAP_TYPE_PCF_SOFT ) || defined( SHADOWMAP_TYPE_VSM )\n\t\t\tvec2 offset = vec2( - 1, 1 ) * shadowRadius * texelSize.y;\n\t\t\treturn (\n\t\t\t\ttexture2DCompare( shadowMap, cubeToUV( bd3D + offset.xyy, texelSize.y ), dp ) +\n\t\t\t\ttexture2DCompare( shadowMap, cubeToUV( bd3D + offset.yyy, texelSize.y ), dp ) +\n\t\t\t\ttexture2DCompare( shadowMap, cubeToUV( bd3D + offset.xyx, texelSize.y ), dp ) +\n\t\t\t\ttexture2DCompare( shadowMap, cubeToUV( bd3D + offset.yyx, texelSize.y ), dp ) +\n\t\t\t\ttexture2DCompare( shadowMap, cubeToUV( bd3D, texelSize.y ), dp ) +\n\t\t\t\ttexture2DCompare( shadowMap, cubeToUV( bd3D + offset.xxy, texelSize.y ), dp ) +\n\t\t\t\ttexture2DCompare( shadowMap, cubeToUV( bd3D + offset.yxy, texelSize.y ), dp ) +\n\t\t\t\ttexture2DCompare( shadowMap, cubeToUV( bd3D + offset.xxx, texelSize.y ), dp ) +\n\t\t\t\ttexture2DCompare( shadowMap, cubeToUV( bd3D + offset.yxx, texelSize.y ), dp )\n\t\t\t) * ( 1.0 / 9.0 );\n\t\t#else\n\t\t\treturn texture2DCompare( shadowMap, cubeToUV( bd3D, texelSize.y ), dp );\n\t\t#endif\n\t}\n#endif",
	shadowmap_pars_vertex = "#ifdef USE_SHADOWMAP\n\t#if NUM_DIR_LIGHT_SHADOWS > 0\n\t\tuniform mat4 directionalShadowMatrix[ NUM_DIR_LIGHT_SHADOWS ];\n\t\tvarying vec4 vDirectionalShadowCoord[ NUM_DIR_LIGHT_SHADOWS ];\n\t\tstruct DirectionalLightShadow {\n\t\t\tfloat shadowBias;\n\t\t\tfloat shadowNormalBias;\n\t\t\tfloat shadowRadius;\n\t\t\tvec2 shadowMapSize;\n\t\t};\n\t\tuniform DirectionalLightShadow directionalLightShadows[ NUM_DIR_LIGHT_SHADOWS ];\n\t#endif\n\t#if NUM_SPOT_LIGHT_SHADOWS > 0\n\t\tuniform mat4 spotShadowMatrix[ NUM_SPOT_LIGHT_SHADOWS ];\n\t\tvarying vec4 vSpotShadowCoord[ NUM_SPOT_LIGHT_SHADOWS ];\n\t\tstruct SpotLightShadow {\n\t\t\tfloat shadowBias;\n\t\t\tfloat shadowNormalBias;\n\t\t\tfloat shadowRadius;\n\t\t\tvec2 shadowMapSize;\n\t\t};\n\t\tuniform SpotLightShadow spotLightShadows[ NUM_SPOT_LIGHT_SHADOWS ];\n\t#endif\n\t#if NUM_POINT_LIGHT_SHADOWS > 0\n\t\tuniform mat4 pointShadowMatrix[ NUM_POINT_LIGHT_SHADOWS ];\n\t\tvarying vec4 vPointShadowCoord[ NUM_POINT_LIGHT_SHADOWS ];\n\t\tstruct PointLightShadow {\n\t\t\tfloat shadowBias;\n\t\t\tfloat shadowNormalBias;\n\t\t\tfloat shadowRadius;\n\t\t\tvec2 shadowMapSize;\n\t\t\tfloat shadowCameraNear;\n\t\t\tfloat shadowCameraFar;\n\t\t};\n\t\tuniform PointLightShadow pointLightShadows[ NUM_POINT_LIGHT_SHADOWS ];\n\t#endif\n#endif",
	shadowmap_vertex = "#ifdef USE_SHADOWMAP\n\t#if NUM_DIR_LIGHT_SHADOWS > 0 || NUM_SPOT_LIGHT_SHADOWS > 0 || NUM_POINT_LIGHT_SHADOWS > 0\n\t\tvec3 shadowWorldNormal = inverseTransformDirection( transformedNormal, viewMatrix );\n\t\tvec4 shadowWorldPosition;\n\t#endif\n\t#if NUM_DIR_LIGHT_SHADOWS > 0\n\t#pragma unroll_loop_start\n\tfor ( int i = 0; i < NUM_DIR_LIGHT_SHADOWS; i ++ ) {\n\t\tshadowWorldPosition = worldPosition + vec4( shadowWorldNormal * directionalLightShadows[ i ].shadowNormalBias, 0 );\n\t\tvDirectionalShadowCoord[ i ] = directionalShadowMatrix[ i ] * shadowWorldPosition;\n\t}\n\t#pragma unroll_loop_end\n\t#endif\n\t#if NUM_SPOT_LIGHT_SHADOWS > 0\n\t#pragma unroll_loop_start\n\tfor ( int i = 0; i < NUM_SPOT_LIGHT_SHADOWS; i ++ ) {\n\t\tshadowWorldPosition = worldPosition + vec4( shadowWorldNormal * spotLightShadows[ i ].shadowNormalBias, 0 );\n\t\tvSpotShadowCoord[ i ] = spotShadowMatrix[ i ] * shadowWorldPosition;\n\t}\n\t#pragma unroll_loop_end\n\t#endif\n\t#if NUM_POINT_LIGHT_SHADOWS > 0\n\t#pragma unroll_loop_start\n\tfor ( int i = 0; i < NUM_POINT_LIGHT_SHADOWS; i ++ ) {\n\t\tshadowWorldPosition = worldPosition + vec4( shadowWorldNormal * pointLightShadows[ i ].shadowNormalBias, 0 );\n\t\tvPointShadowCoord[ i ] = pointShadowMatrix[ i ] * shadowWorldPosition;\n\t}\n\t#pragma unroll_loop_end\n\t#endif\n#endif",
	shadowmask_pars_fragment = "float getShadowMask() {\n\tfloat shadow = 1.0;\n\t#ifdef USE_SHADOWMAP\n\t#if NUM_DIR_LIGHT_SHADOWS > 0\n\tDirectionalLightShadow directionalLight;\n\t#pragma unroll_loop_start\n\tfor ( int i = 0; i < NUM_DIR_LIGHT_SHADOWS; i ++ ) {\n\t\tdirectionalLight = directionalLightShadows[ i ];\n\t\tshadow *= receiveShadow ? getShadow( directionalShadowMap[ i ], directionalLight.shadowMapSize, directionalLight.shadowBias, directionalLight.shadowRadius, vDirectionalShadowCoord[ i ] ) : 1.0;\n\t}\n\t#pragma unroll_loop_end\n\t#endif\n\t#if NUM_SPOT_LIGHT_SHADOWS > 0\n\tSpotLightShadow spotLight;\n\t#pragma unroll_loop_start\n\tfor ( int i = 0; i < NUM_SPOT_LIGHT_SHADOWS; i ++ ) {\n\t\tspotLight = spotLightShadows[ i ];\n\t\tshadow *= receiveShadow ? getShadow( spotShadowMap[ i ], spotLight.shadowMapSize, spotLight.shadowBias, spotLight.shadowRadius, vSpotShadowCoord[ i ] ) : 1.0;\n\t}\n\t#pragma unroll_loop_end\n\t#endif\n\t#if NUM_POINT_LIGHT_SHADOWS > 0\n\tPointLightShadow pointLight;\n\t#pragma unroll_loop_start\n\tfor ( int i = 0; i < NUM_POINT_LIGHT_SHADOWS; i ++ ) {\n\t\tpointLight = pointLightShadows[ i ];\n\t\tshadow *= receiveShadow ? getPointShadow( pointShadowMap[ i ], pointLight.shadowMapSize, pointLight.shadowBias, pointLight.shadowRadius, vPointShadowCoord[ i ], pointLight.shadowCameraNear, pointLight.shadowCameraFar ) : 1.0;\n\t}\n\t#pragma unroll_loop_end\n\t#endif\n\t#endif\n\treturn shadow;\n}",
	skinbase_vertex = "#ifdef USE_SKINNING\n\tmat4 boneMatX = getBoneMatrix( skinIndex.x );\n\tmat4 boneMatY = getBoneMatrix( skinIndex.y );\n\tmat4 boneMatZ = getBoneMatrix( skinIndex.z );\n\tmat4 boneMatW = getBoneMatrix( skinIndex.w );\n#endif",
	skinning_pars_vertex = "#ifdef USE_SKINNING\n\tuniform mat4 bindMatrix;\n\tuniform mat4 bindMatrixInverse;\n\t#ifdef BONE_TEXTURE\n\t\tuniform highp sampler2D boneTexture;\n\t\tuniform int boneTextureSize;\n\t\tmat4 getBoneMatrix( const in float i ) {\n\t\t\tfloat j = i * 4.0;\n\t\t\tfloat x = mod( j, float( boneTextureSize ) );\n\t\t\tfloat y = floor( j / float( boneTextureSize ) );\n\t\t\tfloat dx = 1.0 / float( boneTextureSize );\n\t\t\tfloat dy = 1.0 / float( boneTextureSize );\n\t\t\ty = dy * ( y + 0.5 );\n\t\t\tvec4 v1 = texture2D( boneTexture, vec2( dx * ( x + 0.5 ), y ) );\n\t\t\tvec4 v2 = texture2D( boneTexture, vec2( dx * ( x + 1.5 ), y ) );\n\t\t\tvec4 v3 = texture2D( boneTexture, vec2( dx * ( x + 2.5 ), y ) );\n\t\t\tvec4 v4 = texture2D( boneTexture, vec2( dx * ( x + 3.5 ), y ) );\n\t\t\tmat4 bone = mat4( v1, v2, v3, v4 );\n\t\t\treturn bone;\n\t\t}\n\t#else\n\t\tuniform mat4 boneMatrices[ MAX_BONES ];\n\t\tmat4 getBoneMatrix( const in float i ) {\n\t\t\tmat4 bone = boneMatrices[ int(i) ];\n\t\t\treturn bone;\n\t\t}\n\t#endif\n#endif",
	skinning_vertex = "#ifdef USE_SKINNING\n\tvec4 skinVertex = bindMatrix * vec4( transformed, 1.0 );\n\tvec4 skinned = vec4( 0.0 );\n\tskinned += boneMatX * skinVertex * skinWeight.x;\n\tskinned += boneMatY * skinVertex * skinWeight.y;\n\tskinned += boneMatZ * skinVertex * skinWeight.z;\n\tskinned += boneMatW * skinVertex * skinWeight.w;\n\ttransformed = ( bindMatrixInverse * skinned ).xyz;\n#endif",
	skinnormal_vertex = "#ifdef USE_SKINNING\n\tmat4 skinMatrix = mat4( 0.0 );\n\tskinMatrix += skinWeight.x * boneMatX;\n\tskinMatrix += skinWeight.y * boneMatY;\n\tskinMatrix += skinWeight.z * boneMatZ;\n\tskinMatrix += skinWeight.w * boneMatW;\n\tskinMatrix = bindMatrixInverse * skinMatrix * bindMatrix;\n\tobjectNormal = vec4( skinMatrix * vec4( objectNormal, 0.0 ) ).xyz;\n\t#ifdef USE_TANGENT\n\t\tobjectTangent = vec4( skinMatrix * vec4( objectTangent, 0.0 ) ).xyz;\n\t#endif\n#endif",
	specularmap_fragment = "float specularStrength;\n#ifdef USE_SPECULARMAP\n\tvec4 texelSpecular = texture2D( specularMap, vUv );\n\tspecularStrength = texelSpecular.r;\n#else\n\tspecularStrength = 1.0;\n#endif",
	specularmap_pars_fragment = "#ifdef USE_SPECULARMAP\n\tuniform sampler2D specularMap;\n#endif",
	tonemapping_fragment = "#if defined( TONE_MAPPING )\n\tgl_FragColor.rgb = toneMapping( gl_FragColor.rgb );\n#endif",
	tonemapping_pars_fragment = "#ifndef saturate\n#define saturate(a) clamp( a, 0.0, 1.0 )\n#endif\nuniform float toneMappingExposure;\nvec3 LinearToneMapping( vec3 color ) {\n\treturn toneMappingExposure * color;\n}\nvec3 ReinhardToneMapping( vec3 color ) {\n\tcolor *= toneMappingExposure;\n\treturn saturate( color / ( vec3( 1.0 ) + color ) );\n}\nvec3 OptimizedCineonToneMapping( vec3 color ) {\n\tcolor *= toneMappingExposure;\n\tcolor = max( vec3( 0.0 ), color - 0.004 );\n\treturn pow( ( color * ( 6.2 * color + 0.5 ) ) / ( color * ( 6.2 * color + 1.7 ) + 0.06 ), vec3( 2.2 ) );\n}\nvec3 RRTAndODTFit( vec3 v ) {\n\tvec3 a = v * ( v + 0.0245786 ) - 0.000090537;\n\tvec3 b = v * ( 0.983729 * v + 0.4329510 ) + 0.238081;\n\treturn a / b;\n}\nvec3 ACESFilmicToneMapping( vec3 color ) {\n\tconst mat3 ACESInputMat = mat3(\n\t\tvec3( 0.59719, 0.07600, 0.02840 ),\t\tvec3( 0.35458, 0.90834, 0.13383 ),\n\t\tvec3( 0.04823, 0.01566, 0.83777 )\n\t);\n\tconst mat3 ACESOutputMat = mat3(\n\t\tvec3(  1.60475, -0.10208, -0.00327 ),\t\tvec3( -0.53108,  1.10813, -0.07276 ),\n\t\tvec3( -0.07367, -0.00605,  1.07602 )\n\t);\n\tcolor *= toneMappingExposure / 0.6;\n\tcolor = ACESInputMat * color;\n\tcolor = RRTAndODTFit( color );\n\tcolor = ACESOutputMat * color;\n\treturn saturate( color );\n}\nvec3 CustomToneMapping( vec3 color ) { return color; }",
	transmissionmap_fragment = "#ifdef USE_TRANSMISSIONMAP\n\ttotalTransmission *= texture2D( transmissionMap, vUv ).r;\n#endif",
	transmissionmap_pars_fragment = "#ifdef USE_TRANSMISSIONMAP\n\tuniform sampler2D transmissionMap;\n#endif",
	uv_pars_fragment = "#if ( defined( USE_UV ) && ! defined( UVS_VERTEX_ONLY ) )\n\tvarying vec2 vUv;\n#endif",
	uv_pars_vertex = "#ifdef USE_UV\n\t#ifdef UVS_VERTEX_ONLY\n\t\tvec2 vUv;\n\t#else\n\t\tvarying vec2 vUv;\n\t#endif\n\tuniform mat3 uvTransform;\n#endif",
	uv_vertex = "#ifdef USE_UV\n\tvUv = ( uvTransform * vec3( uv, 1 ) ).xy;\n#endif",
	uv2_pars_fragment = "#if defined( USE_LIGHTMAP ) || defined( USE_AOMAP )\n\tvarying vec2 vUv2;\n#endif",
	uv2_pars_vertex = "#if defined( USE_LIGHTMAP ) || defined( USE_AOMAP )\n\tattribute vec2 uv2;\n\tvarying vec2 vUv2;\n\tuniform mat3 uv2Transform;\n#endif",
	uv2_vertex = "#if defined( USE_LIGHTMAP ) || defined( USE_AOMAP )\n\tvUv2 = ( uv2Transform * vec3( uv2, 1 ) ).xy;\n#endif",
	worldpos_vertex = "#if defined( USE_ENVMAP ) || defined( DISTANCE ) || defined ( USE_SHADOWMAP )\n\tvec4 worldPosition = vec4( transformed, 1.0 );\n\t#ifdef USE_INSTANCING\n\t\tworldPosition = instanceMatrix * worldPosition;\n\t#endif\n\tworldPosition = modelMatrix * worldPosition;\n#endif",
	background_frag = "uniform sampler2D t2D;\nvarying vec2 vUv;\nvoid main() {\n\tvec4 texColor = texture2D( t2D, vUv );\n\tgl_FragColor = mapTexelToLinear( texColor );\n\t#include <tonemapping_fragment>\n\t#include <encodings_fragment>\n}",
	background_vert = "varying vec2 vUv;\nuniform mat3 uvTransform;\nvoid main() {\n\tvUv = ( uvTransform * vec3( uv, 1 ) ).xy;\n\tgl_Position = vec4( position.xy, 1.0, 1.0 );\n}",
	cube_frag = "#include <envmap_common_pars_fragment>\nuniform float opacity;\nvarying vec3 vWorldDirection;\n#include <cube_uv_reflection_fragment>\nvoid main() {\n\tvec3 vReflect = vWorldDirection;\n\t#include <envmap_fragment>\n\tgl_FragColor = envColor;\n\tgl_FragColor.a *= opacity;\n\t#include <tonemapping_fragment>\n\t#include <encodings_fragment>\n}",
	cube_vert = "varying vec3 vWorldDirection;\n#include <common>\nvoid main() {\n\tvWorldDirection = transformDirection( position, modelMatrix );\n\t#include <begin_vertex>\n\t#include <project_vertex>\n\tgl_Position.z = gl_Position.w;\n}",
	depth_frag = "#if DEPTH_PACKING == 3200\n\tuniform float opacity;\n#endif\n#include <common>\n#include <packing>\n#include <uv_pars_fragment>\n#include <map_pars_fragment>\n#include <alphamap_pars_fragment>\n#include <logdepthbuf_pars_fragment>\n#include <clipping_planes_pars_fragment>\nvarying vec2 vHighPrecisionZW;\nvoid main() {\n\t#include <clipping_planes_fragment>\n\tvec4 diffuseColor = vec4( 1.0 );\n\t#if DEPTH_PACKING == 3200\n\t\tdiffuseColor.a = opacity;\n\t#endif\n\t#include <map_fragment>\n\t#include <alphamap_fragment>\n\t#include <alphatest_fragment>\n\t#include <logdepthbuf_fragment>\n\tfloat fragCoordZ = 0.5 * vHighPrecisionZW[0] / vHighPrecisionZW[1] + 0.5;\n\t#if DEPTH_PACKING == 3200\n\t\tgl_FragColor = vec4( vec3( 1.0 - fragCoordZ ), opacity );\n\t#elif DEPTH_PACKING == 3201\n\t\tgl_FragColor = packDepthToRGBA( fragCoordZ );\n\t#endif\n}",
	depth_vert = "#include <common>\n#include <uv_pars_vertex>\n#include <displacementmap_pars_vertex>\n#include <morphtarget_pars_vertex>\n#include <skinning_pars_vertex>\n#include <logdepthbuf_pars_vertex>\n#include <clipping_planes_pars_vertex>\nvarying vec2 vHighPrecisionZW;\nvoid main() {\n\t#include <uv_vertex>\n\t#include <skinbase_vertex>\n\t#ifdef USE_DISPLACEMENTMAP\n\t\t#include <beginnormal_vertex>\n\t\t#include <morphnormal_vertex>\n\t\t#include <skinnormal_vertex>\n\t#endif\n\t#include <begin_vertex>\n\t#include <morphtarget_vertex>\n\t#include <skinning_vertex>\n\t#include <displacementmap_vertex>\n\t#include <project_vertex>\n\t#include <logdepthbuf_vertex>\n\t#include <clipping_planes_vertex>\n\tvHighPrecisionZW = gl_Position.zw;\n}",
	distanceRGBA_frag = "#define DISTANCE\nuniform vec3 referencePosition;\nuniform float nearDistance;\nuniform float farDistance;\nvarying vec3 vWorldPosition;\n#include <common>\n#include <packing>\n#include <uv_pars_fragment>\n#include <map_pars_fragment>\n#include <alphamap_pars_fragment>\n#include <clipping_planes_pars_fragment>\nvoid main () {\n\t#include <clipping_planes_fragment>\n\tvec4 diffuseColor = vec4( 1.0 );\n\t#include <map_fragment>\n\t#include <alphamap_fragment>\n\t#include <alphatest_fragment>\n\tfloat dist = length( vWorldPosition - referencePosition );\n\tdist = ( dist - nearDistance ) / ( farDistance - nearDistance );\n\tdist = saturate( dist );\n\tgl_FragColor = packDepthToRGBA( dist );\n}",
	distanceRGBA_vert = "#define DISTANCE\nvarying vec3 vWorldPosition;\n#include <common>\n#include <uv_pars_vertex>\n#include <displacementmap_pars_vertex>\n#include <morphtarget_pars_vertex>\n#include <skinning_pars_vertex>\n#include <clipping_planes_pars_vertex>\nvoid main() {\n\t#include <uv_vertex>\n\t#include <skinbase_vertex>\n\t#ifdef USE_DISPLACEMENTMAP\n\t\t#include <beginnormal_vertex>\n\t\t#include <morphnormal_vertex>\n\t\t#include <skinnormal_vertex>\n\t#endif\n\t#include <begin_vertex>\n\t#include <morphtarget_vertex>\n\t#include <skinning_vertex>\n\t#include <displacementmap_vertex>\n\t#include <project_vertex>\n\t#include <worldpos_vertex>\n\t#include <clipping_planes_vertex>\n\tvWorldPosition = worldPosition.xyz;\n}",
	equirect_frag = "uniform sampler2D tEquirect;\nvarying vec3 vWorldDirection;\n#include <common>\nvoid main() {\n\tvec3 direction = normalize( vWorldDirection );\n\tvec2 sampleUV = equirectUv( direction );\n\tvec4 texColor = texture2D( tEquirect, sampleUV );\n\tgl_FragColor = mapTexelToLinear( texColor );\n\t#include <tonemapping_fragment>\n\t#include <encodings_fragment>\n}",
	equirect_vert = "varying vec3 vWorldDirection;\n#include <common>\nvoid main() {\n\tvWorldDirection = transformDirection( position, modelMatrix );\n\t#include <begin_vertex>\n\t#include <project_vertex>\n}",
	linedashed_frag = "uniform vec3 diffuse;\nuniform float opacity;\nuniform float dashSize;\nuniform float totalSize;\nvarying float vLineDistance;\n#include <common>\n#include <color_pars_fragment>\n#include <fog_pars_fragment>\n#include <logdepthbuf_pars_fragment>\n#include <clipping_planes_pars_fragment>\nvoid main() {\n\t#include <clipping_planes_fragment>\n\tif ( mod( vLineDistance, totalSize ) > dashSize ) {\n\t\tdiscard;\n\t}\n\tvec3 outgoingLight = vec3( 0.0 );\n\tvec4 diffuseColor = vec4( diffuse, opacity );\n\t#include <logdepthbuf_fragment>\n\t#include <color_fragment>\n\toutgoingLight = diffuseColor.rgb;\n\tgl_FragColor = vec4( outgoingLight, diffuseColor.a );\n\t#include <tonemapping_fragment>\n\t#include <encodings_fragment>\n\t#include <fog_fragment>\n\t#include <premultiplied_alpha_fragment>\n}",
	linedashed_vert = "uniform float scale;\nattribute float lineDistance;\nvarying float vLineDistance;\n#include <common>\n#include <color_pars_vertex>\n#include <fog_pars_vertex>\n#include <morphtarget_pars_vertex>\n#include <logdepthbuf_pars_vertex>\n#include <clipping_planes_pars_vertex>\nvoid main() {\n\tvLineDistance = scale * lineDistance;\n\t#include <color_vertex>\n\t#include <begin_vertex>\n\t#include <morphtarget_vertex>\n\t#include <project_vertex>\n\t#include <logdepthbuf_vertex>\n\t#include <clipping_planes_vertex>\n\t#include <fog_vertex>\n}",
	meshbasic_frag = "uniform vec3 diffuse;\nuniform float opacity;\n#ifndef FLAT_SHADED\n\tvarying vec3 vNormal;\n#endif\n#include <common>\n#include <dithering_pars_fragment>\n#include <color_pars_fragment>\n#include <uv_pars_fragment>\n#include <uv2_pars_fragment>\n#include <map_pars_fragment>\n#include <alphamap_pars_fragment>\n#include <aomap_pars_fragment>\n#include <lightmap_pars_fragment>\n#include <envmap_common_pars_fragment>\n#include <envmap_pars_fragment>\n#include <cube_uv_reflection_fragment>\n#include <fog_pars_fragment>\n#include <specularmap_pars_fragment>\n#include <logdepthbuf_pars_fragment>\n#include <clipping_planes_pars_fragment>\nvoid main() {\n\t#include <clipping_planes_fragment>\n\tvec4 diffuseColor = vec4( diffuse, opacity );\n\t#include <logdepthbuf_fragment>\n\t#include <map_fragment>\n\t#include <color_fragment>\n\t#include <alphamap_fragment>\n\t#include <alphatest_fragment>\n\t#include <specularmap_fragment>\n\tReflectedLight reflectedLight = ReflectedLight( vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ) );\n\t#ifdef USE_LIGHTMAP\n\t\n\t\tvec4 lightMapTexel= texture2D( lightMap, vUv2 );\n\t\treflectedLight.indirectDiffuse += lightMapTexelToLinear( lightMapTexel ).rgb * lightMapIntensity;\n\t#else\n\t\treflectedLight.indirectDiffuse += vec3( 1.0 );\n\t#endif\n\t#include <aomap_fragment>\n\treflectedLight.indirectDiffuse *= diffuseColor.rgb;\n\tvec3 outgoingLight = reflectedLight.indirectDiffuse;\n\t#include <envmap_fragment>\n\tgl_FragColor = vec4( outgoingLight, diffuseColor.a );\n\t#include <tonemapping_fragment>\n\t#include <encodings_fragment>\n\t#include <fog_fragment>\n\t#include <premultiplied_alpha_fragment>\n\t#include <dithering_fragment>\n}",
	meshbasic_vert = "#include <common>\n#include <uv_pars_vertex>\n#include <uv2_pars_vertex>\n#include <envmap_pars_vertex>\n#include <color_pars_vertex>\n#include <fog_pars_vertex>\n#include <morphtarget_pars_vertex>\n#include <skinning_pars_vertex>\n#include <logdepthbuf_pars_vertex>\n#include <clipping_planes_pars_vertex>\nvoid main() {\n\t#include <uv_vertex>\n\t#include <uv2_vertex>\n\t#include <color_vertex>\n\t#include <skinbase_vertex>\n\t#ifdef USE_ENVMAP\n\t#include <beginnormal_vertex>\n\t#include <morphnormal_vertex>\n\t#include <skinnormal_vertex>\n\t#include <defaultnormal_vertex>\n\t#endif\n\t#include <begin_vertex>\n\t#include <morphtarget_vertex>\n\t#include <skinning_vertex>\n\t#include <project_vertex>\n\t#include <logdepthbuf_vertex>\n\t#include <worldpos_vertex>\n\t#include <clipping_planes_vertex>\n\t#include <envmap_vertex>\n\t#include <fog_vertex>\n}",
	meshlambert_frag = "uniform vec3 diffuse;\nuniform vec3 emissive;\nuniform float opacity;\nvarying vec3 vLightFront;\nvarying vec3 vIndirectFront;\n#ifdef DOUBLE_SIDED\n\tvarying vec3 vLightBack;\n\tvarying vec3 vIndirectBack;\n#endif\n#include <common>\n#include <packing>\n#include <dithering_pars_fragment>\n#include <color_pars_fragment>\n#include <uv_pars_fragment>\n#include <uv2_pars_fragment>\n#include <map_pars_fragment>\n#include <alphamap_pars_fragment>\n#include <aomap_pars_fragment>\n#include <lightmap_pars_fragment>\n#include <emissivemap_pars_fragment>\n#include <envmap_common_pars_fragment>\n#include <envmap_pars_fragment>\n#include <cube_uv_reflection_fragment>\n#include <bsdfs>\n#include <lights_pars_begin>\n#include <fog_pars_fragment>\n#include <shadowmap_pars_fragment>\n#include <shadowmask_pars_fragment>\n#include <specularmap_pars_fragment>\n#include <logdepthbuf_pars_fragment>\n#include <clipping_planes_pars_fragment>\nvoid main() {\n\t#include <clipping_planes_fragment>\n\tvec4 diffuseColor = vec4( diffuse, opacity );\n\tReflectedLight reflectedLight = ReflectedLight( vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ) );\n\tvec3 totalEmissiveRadiance = emissive;\n\t#include <logdepthbuf_fragment>\n\t#include <map_fragment>\n\t#include <color_fragment>\n\t#include <alphamap_fragment>\n\t#include <alphatest_fragment>\n\t#include <specularmap_fragment>\n\t#include <emissivemap_fragment>\n\t#ifdef DOUBLE_SIDED\n\t\treflectedLight.indirectDiffuse += ( gl_FrontFacing ) ? vIndirectFront : vIndirectBack;\n\t#else\n\t\treflectedLight.indirectDiffuse += vIndirectFront;\n\t#endif\n\t#include <lightmap_fragment>\n\treflectedLight.indirectDiffuse *= BRDF_Diffuse_Lambert( diffuseColor.rgb );\n\t#ifdef DOUBLE_SIDED\n\t\treflectedLight.directDiffuse = ( gl_FrontFacing ) ? vLightFront : vLightBack;\n\t#else\n\t\treflectedLight.directDiffuse = vLightFront;\n\t#endif\n\treflectedLight.directDiffuse *= BRDF_Diffuse_Lambert( diffuseColor.rgb ) * getShadowMask();\n\t#include <aomap_fragment>\n\tvec3 outgoingLight = reflectedLight.directDiffuse + reflectedLight.indirectDiffuse + totalEmissiveRadiance;\n\t#include <envmap_fragment>\n\tgl_FragColor = vec4( outgoingLight, diffuseColor.a );\n\t#include <tonemapping_fragment>\n\t#include <encodings_fragment>\n\t#include <fog_fragment>\n\t#include <premultiplied_alpha_fragment>\n\t#include <dithering_fragment>\n}",
	meshlambert_vert = "#define LAMBERT\nvarying vec3 vLightFront;\nvarying vec3 vIndirectFront;\n#ifdef DOUBLE_SIDED\n\tvarying vec3 vLightBack;\n\tvarying vec3 vIndirectBack;\n#endif\n#include <common>\n#include <uv_pars_vertex>\n#include <uv2_pars_vertex>\n#include <envmap_pars_vertex>\n#include <bsdfs>\n#include <lights_pars_begin>\n#include <color_pars_vertex>\n#include <fog_pars_vertex>\n#include <morphtarget_pars_vertex>\n#include <skinning_pars_vertex>\n#include <shadowmap_pars_vertex>\n#include <logdepthbuf_pars_vertex>\n#include <clipping_planes_pars_vertex>\nvoid main() {\n\t#include <uv_vertex>\n\t#include <uv2_vertex>\n\t#include <color_vertex>\n\t#include <beginnormal_vertex>\n\t#include <morphnormal_vertex>\n\t#include <skinbase_vertex>\n\t#include <skinnormal_vertex>\n\t#include <defaultnormal_vertex>\n\t#include <begin_vertex>\n\t#include <morphtarget_vertex>\n\t#include <skinning_vertex>\n\t#include <project_vertex>\n\t#include <logdepthbuf_vertex>\n\t#include <clipping_planes_vertex>\n\t#include <worldpos_vertex>\n\t#include <envmap_vertex>\n\t#include <lights_lambert_vertex>\n\t#include <shadowmap_vertex>\n\t#include <fog_vertex>\n}",
	meshmatcap_frag = "#define MATCAP\nuniform vec3 diffuse;\nuniform float opacity;\nuniform sampler2D matcap;\nvarying vec3 vViewPosition;\n#ifndef FLAT_SHADED\n\tvarying vec3 vNormal;\n#endif\n#include <common>\n#include <dithering_pars_fragment>\n#include <color_pars_fragment>\n#include <uv_pars_fragment>\n#include <map_pars_fragment>\n#include <alphamap_pars_fragment>\n#include <fog_pars_fragment>\n#include <bumpmap_pars_fragment>\n#include <normalmap_pars_fragment>\n#include <logdepthbuf_pars_fragment>\n#include <clipping_planes_pars_fragment>\nvoid main() {\n\t#include <clipping_planes_fragment>\n\tvec4 diffuseColor = vec4( diffuse, opacity );\n\t#include <logdepthbuf_fragment>\n\t#include <map_fragment>\n\t#include <color_fragment>\n\t#include <alphamap_fragment>\n\t#include <alphatest_fragment>\n\t#include <normal_fragment_begin>\n\t#include <normal_fragment_maps>\n\tvec3 viewDir = normalize( vViewPosition );\n\tvec3 x = normalize( vec3( viewDir.z, 0.0, - viewDir.x ) );\n\tvec3 y = cross( viewDir, x );\n\tvec2 uv = vec2( dot( x, normal ), dot( y, normal ) ) * 0.495 + 0.5;\n\t#ifdef USE_MATCAP\n\t\tvec4 matcapColor = texture2D( matcap, uv );\n\t\tmatcapColor = matcapTexelToLinear( matcapColor );\n\t#else\n\t\tvec4 matcapColor = vec4( 1.0 );\n\t#endif\n\tvec3 outgoingLight = diffuseColor.rgb * matcapColor.rgb;\n\tgl_FragColor = vec4( outgoingLight, diffuseColor.a );\n\t#include <tonemapping_fragment>\n\t#include <encodings_fragment>\n\t#include <fog_fragment>\n\t#include <premultiplied_alpha_fragment>\n\t#include <dithering_fragment>\n}",
	meshmatcap_vert = "#define MATCAP\nvarying vec3 vViewPosition;\n#ifndef FLAT_SHADED\n\tvarying vec3 vNormal;\n#endif\n#include <common>\n#include <uv_pars_vertex>\n#include <color_pars_vertex>\n#include <displacementmap_pars_vertex>\n#include <fog_pars_vertex>\n#include <morphtarget_pars_vertex>\n#include <skinning_pars_vertex>\n#include <logdepthbuf_pars_vertex>\n#include <clipping_planes_pars_vertex>\nvoid main() {\n\t#include <uv_vertex>\n\t#include <color_vertex>\n\t#include <beginnormal_vertex>\n\t#include <morphnormal_vertex>\n\t#include <skinbase_vertex>\n\t#include <skinnormal_vertex>\n\t#include <defaultnormal_vertex>\n\t#ifndef FLAT_SHADED\n\t\tvNormal = normalize( transformedNormal );\n\t#endif\n\t#include <begin_vertex>\n\t#include <morphtarget_vertex>\n\t#include <skinning_vertex>\n\t#include <displacementmap_vertex>\n\t#include <project_vertex>\n\t#include <logdepthbuf_vertex>\n\t#include <clipping_planes_vertex>\n\t#include <fog_vertex>\n\tvViewPosition = - mvPosition.xyz;\n}",
	meshtoon_frag = "#define TOON\nuniform vec3 diffuse;\nuniform vec3 emissive;\nuniform float opacity;\n#include <common>\n#include <packing>\n#include <dithering_pars_fragment>\n#include <color_pars_fragment>\n#include <uv_pars_fragment>\n#include <uv2_pars_fragment>\n#include <map_pars_fragment>\n#include <alphamap_pars_fragment>\n#include <aomap_pars_fragment>\n#include <lightmap_pars_fragment>\n#include <emissivemap_pars_fragment>\n#include <gradientmap_pars_fragment>\n#include <fog_pars_fragment>\n#include <bsdfs>\n#include <lights_pars_begin>\n#include <lights_toon_pars_fragment>\n#include <shadowmap_pars_fragment>\n#include <bumpmap_pars_fragment>\n#include <normalmap_pars_fragment>\n#include <logdepthbuf_pars_fragment>\n#include <clipping_planes_pars_fragment>\nvoid main() {\n\t#include <clipping_planes_fragment>\n\tvec4 diffuseColor = vec4( diffuse, opacity );\n\tReflectedLight reflectedLight = ReflectedLight( vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ) );\n\tvec3 totalEmissiveRadiance = emissive;\n\t#include <logdepthbuf_fragment>\n\t#include <map_fragment>\n\t#include <color_fragment>\n\t#include <alphamap_fragment>\n\t#include <alphatest_fragment>\n\t#include <normal_fragment_begin>\n\t#include <normal_fragment_maps>\n\t#include <emissivemap_fragment>\n\t#include <lights_toon_fragment>\n\t#include <lights_fragment_begin>\n\t#include <lights_fragment_maps>\n\t#include <lights_fragment_end>\n\t#include <aomap_fragment>\n\tvec3 outgoingLight = reflectedLight.directDiffuse + reflectedLight.indirectDiffuse + totalEmissiveRadiance;\n\tgl_FragColor = vec4( outgoingLight, diffuseColor.a );\n\t#include <tonemapping_fragment>\n\t#include <encodings_fragment>\n\t#include <fog_fragment>\n\t#include <premultiplied_alpha_fragment>\n\t#include <dithering_fragment>\n}",
	meshtoon_vert = "#define TOON\nvarying vec3 vViewPosition;\n#ifndef FLAT_SHADED\n\tvarying vec3 vNormal;\n#endif\n#include <common>\n#include <uv_pars_vertex>\n#include <uv2_pars_vertex>\n#include <displacementmap_pars_vertex>\n#include <color_pars_vertex>\n#include <fog_pars_vertex>\n#include <morphtarget_pars_vertex>\n#include <skinning_pars_vertex>\n#include <shadowmap_pars_vertex>\n#include <logdepthbuf_pars_vertex>\n#include <clipping_planes_pars_vertex>\nvoid main() {\n\t#include <uv_vertex>\n\t#include <uv2_vertex>\n\t#include <color_vertex>\n\t#include <beginnormal_vertex>\n\t#include <morphnormal_vertex>\n\t#include <skinbase_vertex>\n\t#include <skinnormal_vertex>\n\t#include <defaultnormal_vertex>\n#ifndef FLAT_SHADED\n\tvNormal = normalize( transformedNormal );\n#endif\n\t#include <begin_vertex>\n\t#include <morphtarget_vertex>\n\t#include <skinning_vertex>\n\t#include <displacementmap_vertex>\n\t#include <project_vertex>\n\t#include <logdepthbuf_vertex>\n\t#include <clipping_planes_vertex>\n\tvViewPosition = - mvPosition.xyz;\n\t#include <worldpos_vertex>\n\t#include <shadowmap_vertex>\n\t#include <fog_vertex>\n}",
	meshphong_frag = "#define PHONG\nuniform vec3 diffuse;\nuniform vec3 emissive;\nuniform vec3 specular;\nuniform float shininess;\nuniform float opacity;\n#include <common>\n#include <packing>\n#include <dithering_pars_fragment>\n#include <color_pars_fragment>\n#include <uv_pars_fragment>\n#include <uv2_pars_fragment>\n#include <map_pars_fragment>\n#include <alphamap_pars_fragment>\n#include <aomap_pars_fragment>\n#include <lightmap_pars_fragment>\n#include <emissivemap_pars_fragment>\n#include <envmap_common_pars_fragment>\n#include <envmap_pars_fragment>\n#include <cube_uv_reflection_fragment>\n#include <fog_pars_fragment>\n#include <bsdfs>\n#include <lights_pars_begin>\n#include <lights_phong_pars_fragment>\n#include <shadowmap_pars_fragment>\n#include <bumpmap_pars_fragment>\n#include <normalmap_pars_fragment>\n#include <specularmap_pars_fragment>\n#include <logdepthbuf_pars_fragment>\n#include <clipping_planes_pars_fragment>\nvoid main() {\n\t#include <clipping_planes_fragment>\n\tvec4 diffuseColor = vec4( diffuse, opacity );\n\tReflectedLight reflectedLight = ReflectedLight( vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ) );\n\tvec3 totalEmissiveRadiance = emissive;\n\t#include <logdepthbuf_fragment>\n\t#include <map_fragment>\n\t#include <color_fragment>\n\t#include <alphamap_fragment>\n\t#include <alphatest_fragment>\n\t#include <specularmap_fragment>\n\t#include <normal_fragment_begin>\n\t#include <normal_fragment_maps>\n\t#include <emissivemap_fragment>\n\t#include <lights_phong_fragment>\n\t#include <lights_fragment_begin>\n\t#include <lights_fragment_maps>\n\t#include <lights_fragment_end>\n\t#include <aomap_fragment>\n\tvec3 outgoingLight = reflectedLight.directDiffuse + reflectedLight.indirectDiffuse + reflectedLight.directSpecular + reflectedLight.indirectSpecular + totalEmissiveRadiance;\n\t#include <envmap_fragment>\n\tgl_FragColor = vec4( outgoingLight, diffuseColor.a );\n\t#include <tonemapping_fragment>\n\t#include <encodings_fragment>\n\t#include <fog_fragment>\n\t#include <premultiplied_alpha_fragment>\n\t#include <dithering_fragment>\n}",
	meshphong_vert = "#define PHONG\nvarying vec3 vViewPosition;\n#ifndef FLAT_SHADED\n\tvarying vec3 vNormal;\n#endif\n#include <common>\n#include <uv_pars_vertex>\n#include <uv2_pars_vertex>\n#include <displacementmap_pars_vertex>\n#include <envmap_pars_vertex>\n#include <color_pars_vertex>\n#include <fog_pars_vertex>\n#include <morphtarget_pars_vertex>\n#include <skinning_pars_vertex>\n#include <shadowmap_pars_vertex>\n#include <logdepthbuf_pars_vertex>\n#include <clipping_planes_pars_vertex>\nvoid main() {\n\t#include <uv_vertex>\n\t#include <uv2_vertex>\n\t#include <color_vertex>\n\t#include <beginnormal_vertex>\n\t#include <morphnormal_vertex>\n\t#include <skinbase_vertex>\n\t#include <skinnormal_vertex>\n\t#include <defaultnormal_vertex>\n#ifndef FLAT_SHADED\n\tvNormal = normalize( transformedNormal );\n#endif\n\t#include <begin_vertex>\n\t#include <morphtarget_vertex>\n\t#include <skinning_vertex>\n\t#include <displacementmap_vertex>\n\t#include <project_vertex>\n\t#include <logdepthbuf_vertex>\n\t#include <clipping_planes_vertex>\n\tvViewPosition = - mvPosition.xyz;\n\t#include <worldpos_vertex>\n\t#include <envmap_vertex>\n\t#include <shadowmap_vertex>\n\t#include <fog_vertex>\n}",
	meshphysical_frag = "#define STANDARD\n#ifdef PHYSICAL\n\t#define REFLECTIVITY\n\t#define CLEARCOAT\n\t#define TRANSMISSION\n#endif\nuniform vec3 diffuse;\nuniform vec3 emissive;\nuniform float roughness;\nuniform float metalness;\nuniform float opacity;\n#ifdef TRANSMISSION\n\tuniform float transmission;\n#endif\n#ifdef REFLECTIVITY\n\tuniform float reflectivity;\n#endif\n#ifdef CLEARCOAT\n\tuniform float clearcoat;\n\tuniform float clearcoatRoughness;\n#endif\n#ifdef USE_SHEEN\n\tuniform vec3 sheen;\n#endif\nvarying vec3 vViewPosition;\n#ifndef FLAT_SHADED\n\tvarying vec3 vNormal;\n\t#ifdef USE_TANGENT\n\t\tvarying vec3 vTangent;\n\t\tvarying vec3 vBitangent;\n\t#endif\n#endif\n#include <common>\n#include <packing>\n#include <dithering_pars_fragment>\n#include <color_pars_fragment>\n#include <uv_pars_fragment>\n#include <uv2_pars_fragment>\n#include <map_pars_fragment>\n#include <alphamap_pars_fragment>\n#include <aomap_pars_fragment>\n#include <lightmap_pars_fragment>\n#include <emissivemap_pars_fragment>\n#include <transmissionmap_pars_fragment>\n#include <bsdfs>\n#include <cube_uv_reflection_fragment>\n#include <envmap_common_pars_fragment>\n#include <envmap_physical_pars_fragment>\n#include <fog_pars_fragment>\n#include <lights_pars_begin>\n#include <lights_physical_pars_fragment>\n#include <shadowmap_pars_fragment>\n#include <bumpmap_pars_fragment>\n#include <normalmap_pars_fragment>\n#include <clearcoat_pars_fragment>\n#include <roughnessmap_pars_fragment>\n#include <metalnessmap_pars_fragment>\n#include <logdepthbuf_pars_fragment>\n#include <clipping_planes_pars_fragment>\nvoid main() {\n\t#include <clipping_planes_fragment>\n\tvec4 diffuseColor = vec4( diffuse, opacity );\n\tReflectedLight reflectedLight = ReflectedLight( vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ) );\n\tvec3 totalEmissiveRadiance = emissive;\n\t#ifdef TRANSMISSION\n\t\tfloat totalTransmission = transmission;\n\t#endif\n\t#include <logdepthbuf_fragment>\n\t#include <map_fragment>\n\t#include <color_fragment>\n\t#include <alphamap_fragment>\n\t#include <alphatest_fragment>\n\t#include <roughnessmap_fragment>\n\t#include <metalnessmap_fragment>\n\t#include <normal_fragment_begin>\n\t#include <normal_fragment_maps>\n\t#include <clearcoat_normal_fragment_begin>\n\t#include <clearcoat_normal_fragment_maps>\n\t#include <emissivemap_fragment>\n\t#include <transmissionmap_fragment>\n\t#include <lights_physical_fragment>\n\t#include <lights_fragment_begin>\n\t#include <lights_fragment_maps>\n\t#include <lights_fragment_end>\n\t#include <aomap_fragment>\n\tvec3 outgoingLight = reflectedLight.directDiffuse + reflectedLight.indirectDiffuse + reflectedLight.directSpecular + reflectedLight.indirectSpecular + totalEmissiveRadiance;\n\t#ifdef TRANSMISSION\n\t\tdiffuseColor.a *= mix( saturate( 1. - totalTransmission + linearToRelativeLuminance( reflectedLight.directSpecular + reflectedLight.indirectSpecular ) ), 1.0, metalness );\n\t#endif\n\tgl_FragColor = vec4( outgoingLight, diffuseColor.a );\n\t#include <tonemapping_fragment>\n\t#include <encodings_fragment>\n\t#include <fog_fragment>\n\t#include <premultiplied_alpha_fragment>\n\t#include <dithering_fragment>\n}",
	meshphysical_vert = "#define STANDARD\nvarying vec3 vViewPosition;\n#ifndef FLAT_SHADED\n\tvarying vec3 vNormal;\n\t#ifdef USE_TANGENT\n\t\tvarying vec3 vTangent;\n\t\tvarying vec3 vBitangent;\n\t#endif\n#endif\n#include <common>\n#include <uv_pars_vertex>\n#include <uv2_pars_vertex>\n#include <displacementmap_pars_vertex>\n#include <color_pars_vertex>\n#include <fog_pars_vertex>\n#include <morphtarget_pars_vertex>\n#include <skinning_pars_vertex>\n#include <shadowmap_pars_vertex>\n#include <logdepthbuf_pars_vertex>\n#include <clipping_planes_pars_vertex>\nvoid main() {\n\t#include <uv_vertex>\n\t#include <uv2_vertex>\n\t#include <color_vertex>\n\t#include <beginnormal_vertex>\n\t#include <morphnormal_vertex>\n\t#include <skinbase_vertex>\n\t#include <skinnormal_vertex>\n\t#include <defaultnormal_vertex>\n#ifndef FLAT_SHADED\n\tvNormal = normalize( transformedNormal );\n\t#ifdef USE_TANGENT\n\t\tvTangent = normalize( transformedTangent );\n\t\tvBitangent = normalize( cross( vNormal, vTangent ) * tangent.w );\n\t#endif\n#endif\n\t#include <begin_vertex>\n\t#include <morphtarget_vertex>\n\t#include <skinning_vertex>\n\t#include <displacementmap_vertex>\n\t#include <project_vertex>\n\t#include <logdepthbuf_vertex>\n\t#include <clipping_planes_vertex>\n\tvViewPosition = - mvPosition.xyz;\n\t#include <worldpos_vertex>\n\t#include <shadowmap_vertex>\n\t#include <fog_vertex>\n}",
	normal_frag = "#define NORMAL\nuniform float opacity;\n#if defined( FLAT_SHADED ) || defined( USE_BUMPMAP ) || defined( TANGENTSPACE_NORMALMAP )\n\tvarying vec3 vViewPosition;\n#endif\n#ifndef FLAT_SHADED\n\tvarying vec3 vNormal;\n\t#ifdef USE_TANGENT\n\t\tvarying vec3 vTangent;\n\t\tvarying vec3 vBitangent;\n\t#endif\n#endif\n#include <packing>\n#include <uv_pars_fragment>\n#include <bumpmap_pars_fragment>\n#include <normalmap_pars_fragment>\n#include <logdepthbuf_pars_fragment>\n#include <clipping_planes_pars_fragment>\nvoid main() {\n\t#include <clipping_planes_fragment>\n\t#include <logdepthbuf_fragment>\n\t#include <normal_fragment_begin>\n\t#include <normal_fragment_maps>\n\tgl_FragColor = vec4( packNormalToRGB( normal ), opacity );\n}",
	normal_vert = "#define NORMAL\n#if defined( FLAT_SHADED ) || defined( USE_BUMPMAP ) || defined( TANGENTSPACE_NORMALMAP )\n\tvarying vec3 vViewPosition;\n#endif\n#ifndef FLAT_SHADED\n\tvarying vec3 vNormal;\n\t#ifdef USE_TANGENT\n\t\tvarying vec3 vTangent;\n\t\tvarying vec3 vBitangent;\n\t#endif\n#endif\n#include <common>\n#include <uv_pars_vertex>\n#include <displacementmap_pars_vertex>\n#include <morphtarget_pars_vertex>\n#include <skinning_pars_vertex>\n#include <logdepthbuf_pars_vertex>\n#include <clipping_planes_pars_vertex>\nvoid main() {\n\t#include <uv_vertex>\n\t#include <beginnormal_vertex>\n\t#include <morphnormal_vertex>\n\t#include <skinbase_vertex>\n\t#include <skinnormal_vertex>\n\t#include <defaultnormal_vertex>\n#ifndef FLAT_SHADED\n\tvNormal = normalize( transformedNormal );\n\t#ifdef USE_TANGENT\n\t\tvTangent = normalize( transformedTangent );\n\t\tvBitangent = normalize( cross( vNormal, vTangent ) * tangent.w );\n\t#endif\n#endif\n\t#include <begin_vertex>\n\t#include <morphtarget_vertex>\n\t#include <skinning_vertex>\n\t#include <displacementmap_vertex>\n\t#include <project_vertex>\n\t#include <logdepthbuf_vertex>\n\t#include <clipping_planes_vertex>\n#if defined( FLAT_SHADED ) || defined( USE_BUMPMAP ) || defined( TANGENTSPACE_NORMALMAP )\n\tvViewPosition = - mvPosition.xyz;\n#endif\n}",
	points_frag = "uniform vec3 diffuse;\nuniform float opacity;\n#include <common>\n#include <color_pars_fragment>\n#include <map_particle_pars_fragment>\n#include <fog_pars_fragment>\n#include <logdepthbuf_pars_fragment>\n#include <clipping_planes_pars_fragment>\nvoid main() {\n\t#include <clipping_planes_fragment>\n\tvec3 outgoingLight = vec3( 0.0 );\n\tvec4 diffuseColor = vec4( diffuse, opacity );\n\t#include <logdepthbuf_fragment>\n\t#include <map_particle_fragment>\n\t#include <color_fragment>\n\t#include <alphatest_fragment>\n\toutgoingLight = diffuseColor.rgb;\n\tgl_FragColor = vec4( outgoingLight, diffuseColor.a );\n\t#include <tonemapping_fragment>\n\t#include <encodings_fragment>\n\t#include <fog_fragment>\n\t#include <premultiplied_alpha_fragment>\n}",
	points_vert = "uniform float size;\nuniform float scale;\n#include <common>\n#include <color_pars_vertex>\n#include <fog_pars_vertex>\n#include <morphtarget_pars_vertex>\n#include <logdepthbuf_pars_vertex>\n#include <clipping_planes_pars_vertex>\nvoid main() {\n\t#include <color_vertex>\n\t#include <begin_vertex>\n\t#include <morphtarget_vertex>\n\t#include <project_vertex>\n\tgl_PointSize = size;\n\t#ifdef USE_SIZEATTENUATION\n\t\tbool isPerspective = isPerspectiveMatrix( projectionMatrix );\n\t\tif ( isPerspective ) gl_PointSize *= ( scale / - mvPosition.z );\n\t#endif\n\t#include <logdepthbuf_vertex>\n\t#include <clipping_planes_vertex>\n\t#include <worldpos_vertex>\n\t#include <fog_vertex>\n}",
	shadow_frag = "uniform vec3 color;\nuniform float opacity;\n#include <common>\n#include <packing>\n#include <fog_pars_fragment>\n#include <bsdfs>\n#include <lights_pars_begin>\n#include <shadowmap_pars_fragment>\n#include <shadowmask_pars_fragment>\nvoid main() {\n\tgl_FragColor = vec4( color, opacity * ( 1.0 - getShadowMask() ) );\n\t#include <tonemapping_fragment>\n\t#include <encodings_fragment>\n\t#include <fog_fragment>\n}",
	shadow_vert = "#include <common>\n#include <fog_pars_vertex>\n#include <shadowmap_pars_vertex>\nvoid main() {\n\t#include <begin_vertex>\n\t#include <project_vertex>\n\t#include <worldpos_vertex>\n\t#include <beginnormal_vertex>\n\t#include <morphnormal_vertex>\n\t#include <skinbase_vertex>\n\t#include <skinnormal_vertex>\n\t#include <defaultnormal_vertex>\n\t#include <shadowmap_vertex>\n\t#include <fog_vertex>\n}",
	sprite_frag = "uniform vec3 diffuse;\nuniform float opacity;\n#include <common>\n#include <uv_pars_fragment>\n#include <map_pars_fragment>\n#include <alphamap_pars_fragment>\n#include <fog_pars_fragment>\n#include <logdepthbuf_pars_fragment>\n#include <clipping_planes_pars_fragment>\nvoid main() {\n\t#include <clipping_planes_fragment>\n\tvec3 outgoingLight = vec3( 0.0 );\n\tvec4 diffuseColor = vec4( diffuse, opacity );\n\t#include <logdepthbuf_fragment>\n\t#include <map_fragment>\n\t#include <alphamap_fragment>\n\t#include <alphatest_fragment>\n\toutgoingLight = diffuseColor.rgb;\n\tgl_FragColor = vec4( outgoingLight, diffuseColor.a );\n\t#include <tonemapping_fragment>\n\t#include <encodings_fragment>\n\t#include <fog_fragment>\n}",
	sprite_vert = "uniform float rotation;\nuniform vec2 center;\n#include <common>\n#include <uv_pars_vertex>\n#include <fog_pars_vertex>\n#include <logdepthbuf_pars_vertex>\n#include <clipping_planes_pars_vertex>\nvoid main() {\n\t#include <uv_vertex>\n\tvec4 mvPosition = modelViewMatrix * vec4( 0.0, 0.0, 0.0, 1.0 );\n\tvec2 scale;\n\tscale.x = length( vec3( modelMatrix[ 0 ].x, modelMatrix[ 0 ].y, modelMatrix[ 0 ].z ) );\n\tscale.y = length( vec3( modelMatrix[ 1 ].x, modelMatrix[ 1 ].y, modelMatrix[ 1 ].z ) );\n\t#ifndef USE_SIZEATTENUATION\n\t\tbool isPerspective = isPerspectiveMatrix( projectionMatrix );\n\t\tif ( isPerspective ) scale *= - mvPosition.z;\n\t#endif\n\tvec2 alignedPosition = ( position.xy - ( center - vec2( 0.5 ) ) ) * scale;\n\tvec2 rotatedPosition;\n\trotatedPosition.x = cos( rotation ) * alignedPosition.x - sin( rotation ) * alignedPosition.y;\n\trotatedPosition.y = sin( rotation ) * alignedPosition.x + cos( rotation ) * alignedPosition.y;\n\tmvPosition.xy += rotatedPosition;\n\tgl_Position = projectionMatrix * mvPosition;\n\t#include <logdepthbuf_vertex>\n\t#include <clipping_planes_vertex>\n\t#include <fog_vertex>\n}";
const ShaderChunk = {
		alphamap_fragment: alphamap_fragment,
		alphamap_pars_fragment: alphamap_pars_fragment,
		alphatest_fragment: alphatest_fragment,
		aomap_fragment: aomap_fragment,
		aomap_pars_fragment: aomap_pars_fragment,
		begin_vertex: begin_vertex,
		beginnormal_vertex: beginnormal_vertex,
		bsdfs: bsdfs,
		bumpmap_pars_fragment: bumpmap_pars_fragment,
		clipping_planes_fragment: clipping_planes_fragment,
		clipping_planes_pars_fragment: clipping_planes_pars_fragment,
		clipping_planes_pars_vertex: clipping_planes_pars_vertex,
		clipping_planes_vertex: clipping_planes_vertex,
		color_fragment: color_fragment,
		color_pars_fragment: color_pars_fragment,
		color_pars_vertex: color_pars_vertex,
		color_vertex: color_vertex,
		common: common,
		cube_uv_reflection_fragment: cube_uv_reflection_fragment,
		defaultnormal_vertex: defaultnormal_vertex,
		displacementmap_pars_vertex: displacementmap_pars_vertex,
		displacementmap_vertex: displacementmap_vertex,
		emissivemap_fragment: emissivemap_fragment,
		emissivemap_pars_fragment: emissivemap_pars_fragment,
		encodings_fragment: encodings_fragment,
		encodings_pars_fragment: encodings_pars_fragment,
		envmap_fragment: envmap_fragment,
		envmap_common_pars_fragment: envmap_common_pars_fragment,
		envmap_pars_fragment: envmap_pars_fragment,
		envmap_pars_vertex: envmap_pars_vertex,
		envmap_physical_pars_fragment: envmap_physical_pars_fragment,
		envmap_vertex: envmap_vertex,
		fog_vertex: fog_vertex,
		fog_pars_vertex: fog_pars_vertex,
		fog_fragment: fog_fragment,
		fog_pars_fragment: fog_pars_fragment,
		gradientmap_pars_fragment: gradientmap_pars_fragment,
		lightmap_fragment: lightmap_fragment,
		lightmap_pars_fragment: lightmap_pars_fragment,
		lights_lambert_vertex: lights_lambert_vertex,
		lights_pars_begin: lights_pars_begin,
		lights_toon_fragment: lights_toon_fragment,
		lights_toon_pars_fragment: lights_toon_pars_fragment,
		lights_phong_fragment: lights_phong_fragment,
		lights_phong_pars_fragment: lights_phong_pars_fragment,
		lights_physical_fragment: lights_physical_fragment,
		lights_physical_pars_fragment: lights_physical_pars_fragment,
		lights_fragment_begin: lights_fragment_begin,
		lights_fragment_maps: lights_fragment_maps,
		lights_fragment_end: lights_fragment_end,
		logdepthbuf_fragment: logdepthbuf_fragment,
		logdepthbuf_pars_fragment: logdepthbuf_pars_fragment,
		logdepthbuf_pars_vertex: logdepthbuf_pars_vertex,
		logdepthbuf_vertex: logdepthbuf_vertex,
		map_fragment: map_fragment,
		map_pars_fragment: map_pars_fragment,
		map_particle_fragment: map_particle_fragment,
		map_particle_pars_fragment: map_particle_pars_fragment,
		metalnessmap_fragment: metalnessmap_fragment,
		metalnessmap_pars_fragment: metalnessmap_pars_fragment,
		morphnormal_vertex: morphnormal_vertex,
		morphtarget_pars_vertex: morphtarget_pars_vertex,
		morphtarget_vertex: morphtarget_vertex,
		normal_fragment_begin: normal_fragment_begin,
		normal_fragment_maps: normal_fragment_maps,
		normalmap_pars_fragment: normalmap_pars_fragment,
		clearcoat_normal_fragment_begin: clearcoat_normal_fragment_begin,
		clearcoat_normal_fragment_maps: clearcoat_normal_fragment_maps,
		clearcoat_pars_fragment: clearcoat_pars_fragment,
		packing: packing,
		premultiplied_alpha_fragment: premultiplied_alpha_fragment,
		project_vertex: project_vertex,
		dithering_fragment: dithering_fragment,
		dithering_pars_fragment: dithering_pars_fragment,
		roughnessmap_fragment: roughnessmap_fragment,
		roughnessmap_pars_fragment: roughnessmap_pars_fragment,
		shadowmap_pars_fragment: shadowmap_pars_fragment,
		shadowmap_pars_vertex: shadowmap_pars_vertex,
		shadowmap_vertex: shadowmap_vertex,
		shadowmask_pars_fragment: shadowmask_pars_fragment,
		skinbase_vertex: skinbase_vertex,
		skinning_pars_vertex: skinning_pars_vertex,
		skinning_vertex: skinning_vertex,
		skinnormal_vertex: skinnormal_vertex,
		specularmap_fragment: specularmap_fragment,
		specularmap_pars_fragment: specularmap_pars_fragment,
		tonemapping_fragment: tonemapping_fragment,
		tonemapping_pars_fragment: tonemapping_pars_fragment,
		transmissionmap_fragment: transmissionmap_fragment,
		transmissionmap_pars_fragment: transmissionmap_pars_fragment,
		uv_pars_fragment: uv_pars_fragment,
		uv_pars_vertex: uv_pars_vertex,
		uv_vertex: uv_vertex,
		uv2_pars_fragment: uv2_pars_fragment,
		uv2_pars_vertex: uv2_pars_vertex,
		uv2_vertex: uv2_vertex,
		worldpos_vertex: worldpos_vertex,
		background_frag: background_frag,
		background_vert: background_vert,
		cube_frag: cube_frag,
		cube_vert: cube_vert,
		depth_frag: depth_frag,
		depth_vert: depth_vert,
		distanceRGBA_frag: distanceRGBA_frag,
		distanceRGBA_vert: distanceRGBA_vert,
		equirect_frag: equirect_frag,
		equirect_vert: equirect_vert,
		linedashed_frag: linedashed_frag,
		linedashed_vert: linedashed_vert,
		meshbasic_frag: meshbasic_frag,
		meshbasic_vert: meshbasic_vert,
		meshlambert_frag: meshlambert_frag,
		meshlambert_vert: meshlambert_vert,
		meshmatcap_frag: meshmatcap_frag,
		meshmatcap_vert: meshmatcap_vert,
		meshtoon_frag: meshtoon_frag,
		meshtoon_vert: meshtoon_vert,
		meshphong_frag: meshphong_frag,
		meshphong_vert: meshphong_vert,
		meshphysical_frag: meshphysical_frag,
		meshphysical_vert: meshphysical_vert,
		normal_frag: normal_frag,
		normal_vert: normal_vert,
		points_frag: points_frag,
		points_vert: points_vert,
		shadow_frag: shadow_frag,
		shadow_vert: shadow_vert,
		sprite_frag: sprite_frag,
		sprite_vert: sprite_vert
	},
	UniformsLib = {
		common: {
			diffuse: {
				value: new Color(15658734)
			},
			opacity: {
				value: 1
			},
			map: {
				value: null
			},
			uvTransform: {
				value: new Matrix3
			},
			uv2Transform: {
				value: new Matrix3
			},
			alphaMap: {
				value: null
			}
		},
		specularmap: {
			specularMap: {
				value: null
			}
		},
		envmap: {
			envMap: {
				value: null
			},
			flipEnvMap: {
				value: -1
			},
			reflectivity: {
				value: 1
			},
			refractionRatio: {
				value: .98
			},
			maxMipLevel: {
				value: 0
			}
		},
		aomap: {
			aoMap: {
				value: null
			},
			aoMapIntensity: {
				value: 1
			}
		},
		lightmap: {
			lightMap: {
				value: null
			},
			lightMapIntensity: {
				value: 1
			}
		},
		emissivemap: {
			emissiveMap: {
				value: null
			}
		},
		bumpmap: {
			bumpMap: {
				value: null
			},
			bumpScale: {
				value: 1
			}
		},
		normalmap: {
			normalMap: {
				value: null
			},
			normalScale: {
				value: new Vector2(1, 1)
			}
		},
		displacementmap: {
			displacementMap: {
				value: null
			},
			displacementScale: {
				value: 1
			},
			displacementBias: {
				value: 0
			}
		},
		roughnessmap: {
			roughnessMap: {
				value: null
			}
		},
		metalnessmap: {
			metalnessMap: {
				value: null
			}
		},
		gradientmap: {
			gradientMap: {
				value: null
			}
		},
		fog: {
			fogDensity: {
				value: 25e-5
			},
			fogNear: {
				value: 1
			},
			fogFar: {
				value: 2e3
			},
			fogColor: {
				value: new Color(16777215)
			}
		},
		lights: {
			ambientLightColor: {
				value: []
			},
			lightProbe: {
				value: []
			},
			directionalLights: {
				value: [],
				properties: {
					direction: {},
					color: {}
				}
			},
			directionalLightShadows: {
				value: [],
				properties: {
					shadowBias: {},
					shadowNormalBias: {},
					shadowRadius: {},
					shadowMapSize: {}
				}
			},
			directionalShadowMap: {
				value: []
			},
			directionalShadowMatrix: {
				value: []
			},
			spotLights: {
				value: [],
				properties: {
					color: {},
					position: {},
					direction: {},
					distance: {},
					coneCos: {},
					penumbraCos: {},
					decay: {}
				}
			},
			spotLightShadows: {
				value: [],
				properties: {
					shadowBias: {},
					shadowNormalBias: {},
					shadowRadius: {},
					shadowMapSize: {}
				}
			},
			spotShadowMap: {
				value: []
			},
			spotShadowMatrix: {
				value: []
			},
			pointLights: {
				value: [],
				properties: {
					color: {},
					position: {},
					decay: {},
					distance: {}
				}
			},
			pointLightShadows: {
				value: [],
				properties: {
					shadowBias: {},
					shadowNormalBias: {},
					shadowRadius: {},
					shadowMapSize: {},
					shadowCameraNear: {},
					shadowCameraFar: {}
				}
			},
			pointShadowMap: {
				value: []
			},
			pointShadowMatrix: {
				value: []
			},
			hemisphereLights: {
				value: [],
				properties: {
					direction: {},
					skyColor: {},
					groundColor: {}
				}
			},
			rectAreaLights: {
				value: [],
				properties: {
					color: {},
					position: {},
					width: {},
					height: {}
				}
			},
			ltc_1: {
				value: null
			},
			ltc_2: {
				value: null
			}
		},
		points: {
			diffuse: {
				value: new Color(15658734)
			},
			opacity: {
				value: 1
			},
			size: {
				value: 1
			},
			scale: {
				value: 1
			},
			map: {
				value: null
			},
			alphaMap: {
				value: null
			},
			uvTransform: {
				value: new Matrix3
			}
		},
		sprite: {
			diffuse: {
				value: new Color(15658734)
			},
			opacity: {
				value: 1
			},
			center: {
				value: new Vector2(.5, .5)
			},
			rotation: {
				value: 0
			},
			map: {
				value: null
			},
			alphaMap: {
				value: null
			},
			uvTransform: {
				value: new Matrix3
			}
		}
	},
	ShaderLib = {
		basic: {
			uniforms: mergeUniforms([UniformsLib.common, UniformsLib.specularmap, UniformsLib.envmap, UniformsLib.aomap, UniformsLib.lightmap, UniformsLib.fog]),
			vertexShader: ShaderChunk.meshbasic_vert,
			fragmentShader: ShaderChunk.meshbasic_frag
		},
		lambert: {
			uniforms: mergeUniforms([UniformsLib.common, UniformsLib.specularmap, UniformsLib.envmap, UniformsLib.aomap, UniformsLib.lightmap, UniformsLib.emissivemap, UniformsLib.fog, UniformsLib.lights, {
				emissive: {
					value: new Color(0)
				}
			}]),
			vertexShader: ShaderChunk.meshlambert_vert,
			fragmentShader: ShaderChunk.meshlambert_frag
		},
		phong: {
			uniforms: mergeUniforms([UniformsLib.common, UniformsLib.specularmap, UniformsLib.envmap, UniformsLib.aomap, UniformsLib.lightmap, UniformsLib.emissivemap, UniformsLib.bumpmap, UniformsLib.normalmap, UniformsLib.displacementmap, UniformsLib.fog, UniformsLib.lights, {
				emissive: {
					value: new Color(0)
				},
				specular: {
					value: new Color(1118481)
				},
				shininess: {
					value: 30
				}
			}]),
			vertexShader: ShaderChunk.meshphong_vert,
			fragmentShader: ShaderChunk.meshphong_frag
		},
		standard: {
			uniforms: mergeUniforms([UniformsLib.common, UniformsLib.envmap, UniformsLib.aomap, UniformsLib.lightmap, UniformsLib.emissivemap, UniformsLib.bumpmap, UniformsLib.normalmap, UniformsLib.displacementmap, UniformsLib.roughnessmap, UniformsLib.metalnessmap, UniformsLib.fog, UniformsLib.lights, {
				emissive: {
					value: new Color(0)
				},
				roughness: {
					value: 1
				},
				metalness: {
					value: 0
				},
				envMapIntensity: {
					value: 1
				}
			}]),
			vertexShader: ShaderChunk.meshphysical_vert,
			fragmentShader: ShaderChunk.meshphysical_frag
		},
		toon: {
			uniforms: mergeUniforms([UniformsLib.common, UniformsLib.aomap, UniformsLib.lightmap, UniformsLib.emissivemap, UniformsLib.bumpmap, UniformsLib.normalmap, UniformsLib.displacementmap, UniformsLib.gradientmap, UniformsLib.fog, UniformsLib.lights, {
				emissive: {
					value: new Color(0)
				}
			}]),
			vertexShader: ShaderChunk.meshtoon_vert,
			fragmentShader: ShaderChunk.meshtoon_frag
		},
		matcap: {
			uniforms: mergeUniforms([UniformsLib.common, UniformsLib.bumpmap, UniformsLib.normalmap, UniformsLib.displacementmap, UniformsLib.fog, {
				matcap: {
					value: null
				}
			}]),
			vertexShader: ShaderChunk.meshmatcap_vert,
			fragmentShader: ShaderChunk.meshmatcap_frag
		},
		points: {
			uniforms: mergeUniforms([UniformsLib.points, UniformsLib.fog]),
			vertexShader: ShaderChunk.points_vert,
			fragmentShader: ShaderChunk.points_frag
		},
		dashed: {
			uniforms: mergeUniforms([UniformsLib.common, UniformsLib.fog, {
				scale: {
					value: 1
				},
				dashSize: {
					value: 1
				},
				totalSize: {
					value: 2
				}
			}]),
			vertexShader: ShaderChunk.linedashed_vert,
			fragmentShader: ShaderChunk.linedashed_frag
		},
		depth: {
			uniforms: mergeUniforms([UniformsLib.common, UniformsLib.displacementmap]),
			vertexShader: ShaderChunk.depth_vert,
			fragmentShader: ShaderChunk.depth_frag
		},
		normal: {
			uniforms: mergeUniforms([UniformsLib.common, UniformsLib.bumpmap, UniformsLib.normalmap, UniformsLib.displacementmap, {
				opacity: {
					value: 1
				}
			}]),
			vertexShader: ShaderChunk.normal_vert,
			fragmentShader: ShaderChunk.normal_frag
		},
		sprite: {
			uniforms: mergeUniforms([UniformsLib.sprite, UniformsLib.fog]),
			vertexShader: ShaderChunk.sprite_vert,
			fragmentShader: ShaderChunk.sprite_frag
		},
		background: {
			uniforms: {
				uvTransform: {
					value: new Matrix3
				},
				t2D: {
					value: null
				}
			},
			vertexShader: ShaderChunk.background_vert,
			fragmentShader: ShaderChunk.background_frag
		},
		cube: {
			uniforms: mergeUniforms([UniformsLib.envmap, {
				opacity: {
					value: 1
				}
			}]),
			vertexShader: ShaderChunk.cube_vert,
			fragmentShader: ShaderChunk.cube_frag
		},
		equirect: {
			uniforms: {
				tEquirect: {
					value: null
				}
			},
			vertexShader: ShaderChunk.equirect_vert,
			fragmentShader: ShaderChunk.equirect_frag
		},
		distanceRGBA: {
			uniforms: mergeUniforms([UniformsLib.common, UniformsLib.displacementmap, {
				referencePosition: {
					value: new Vector3
				},
				nearDistance: {
					value: 1
				},
				farDistance: {
					value: 1e3
				}
			}]),
			vertexShader: ShaderChunk.distanceRGBA_vert,
			fragmentShader: ShaderChunk.distanceRGBA_frag
		},
		shadow: {
			uniforms: mergeUniforms([UniformsLib.lights, UniformsLib.fog, {
				color: {
					value: new Color(0)
				},
				opacity: {
					value: 1
				}
			}]),
			vertexShader: ShaderChunk.shadow_vert,
			fragmentShader: ShaderChunk.shadow_frag
		}
	};

function WebGLBackground(e, t, i, n, r) {
	const a = new Color(0);
	let s, o, l = 0,
		c = null,
		h = 0,
		u = null;

	function d(e, t) {
		i.buffers.color.setClear(e.r, e.g, e.b, t, r)
	}
	return {
		getClearColor: function() {
			return a
		},
		setClearColor: function(e, t = 1) {
			a.set(e), l = t, d(a, l)
		},
		getClearAlpha: function() {
			return l
		},
		setClearAlpha: function(e) {
			l = e, d(a, l)
		},
		render: function(i, r, p, m) {
			let A = !0 === r.isScene ? r.background : null;
			A && A.isTexture && (A = t.get(A));
			const g = e.xr,
				f = g.getSession && g.getSession();
			f && "additive" === f.environmentBlendMode && (A = null), null === A ? d(a, l) : A && A.isColor && (d(A, 1), m = !0), (e.autoClear || m) && e.clear(e.autoClearColor, e.autoClearDepth, e.autoClearStencil), A && (A.isCubeTexture || A.isWebGLCubeRenderTarget || 306 === A.mapping) ? (void 0 === o && (o = new Mesh(new BoxGeometry(1, 1, 1), new ShaderMaterial({
				name: "BackgroundCubeMaterial",
				uniforms: cloneUniforms(ShaderLib.cube.uniforms),
				vertexShader: ShaderLib.cube.vertexShader,
				fragmentShader: ShaderLib.cube.fragmentShader,
				side: 1,
				depthTest: !1,
				depthWrite: !1,
				fog: !1
			})), o.geometry.deleteAttribute("normal"), o.geometry.deleteAttribute("uv"), o.onBeforeRender = function(e, t, i) {
				this.matrixWorld.copyPosition(i.matrixWorld)
			}, Object.defineProperty(o.material, "envMap", {
				get: function() {
					return this.uniforms.envMap.value
				}
			}), n.update(o)), A.isWebGLCubeRenderTarget && (A = A.texture), o.material.uniforms.envMap.value = A, o.material.uniforms.flipEnvMap.value = A.isCubeTexture && A._needsFlipEnvMap ? -1 : 1, c === A && h === A.version && u === e.toneMapping || (o.material.needsUpdate = !0, c = A, h = A.version, u = e.toneMapping), i.unshift(o, o.geometry, o.material, 0, 0, null)) : A && A.isTexture && (void 0 === s && (s = new Mesh(new PlaneGeometry(2, 2), new ShaderMaterial({
				name: "BackgroundMaterial",
				uniforms: cloneUniforms(ShaderLib.background.uniforms),
				vertexShader: ShaderLib.background.vertexShader,
				fragmentShader: ShaderLib.background.fragmentShader,
				side: 0,
				depthTest: !1,
				depthWrite: !1,
				fog: !1
			})), s.geometry.deleteAttribute("normal"), Object.defineProperty(s.material, "map", {
				get: function() {
					return this.uniforms.t2D.value
				}
			}), n.update(s)), s.material.uniforms.t2D.value = A, !0 === A.matrixAutoUpdate && A.updateMatrix(), s.material.uniforms.uvTransform.value.copy(A.matrix), c === A && h === A.version && u === e.toneMapping || (s.material.needsUpdate = !0, c = A, h = A.version, u = e.toneMapping), i.unshift(s, s.geometry, s.material, 0, 0, null))
		}
	}
}

function WebGLBindingStates(e, t, i, n) {
	const r = e.getParameter(34921),
		a = n.isWebGL2 ? null : t.get("OES_vertex_array_object"),
		s = n.isWebGL2 || null !== a,
		o = {},
		l = d(null);
	let c = l;

	function h(t) {
		return n.isWebGL2 ? e.bindVertexArray(t) : a.bindVertexArrayOES(t)
	}

	function u(t) {
		return n.isWebGL2 ? e.deleteVertexArray(t) : a.deleteVertexArrayOES(t)
	}

	function d(e) {
		const t = [],
			i = [],
			n = [];
		for (let e = 0; e < r; e++) t[e] = 0, i[e] = 0, n[e] = 0;
		return {
			geometry: null,
			program: null,
			wireframe: !1,
			newAttributes: t,
			enabledAttributes: i,
			attributeDivisors: n,
			object: e,
			attributes: {},
			index: null
		}
	}

	function p() {
		const e = c.newAttributes;
		for (let t = 0, i = e.length; t < i; t++) e[t] = 0
	}

	function m(e) {
		A(e, 0)
	}

	function A(i, r) {
		const a = c.newAttributes,
			s = c.enabledAttributes,
			o = c.attributeDivisors;
		if (a[i] = 1, 0 === s[i] && (e.enableVertexAttribArray(i), s[i] = 1), o[i] !== r) {
			(n.isWebGL2 ? e : t.get("ANGLE_instanced_arrays"))[n.isWebGL2 ? "vertexAttribDivisor" : "vertexAttribDivisorANGLE"](i, r), o[i] = r
		}
	}

	function g() {
		const t = c.newAttributes,
			i = c.enabledAttributes;
		for (let n = 0, r = i.length; n < r; n++) i[n] !== t[n] && (e.disableVertexAttribArray(n), i[n] = 0)
	}

	function f(t, i, r, a, s, o) {
		!0 !== n.isWebGL2 || 5124 !== r && 5125 !== r ? e.vertexAttribPointer(t, i, r, a, s, o) : e.vertexAttribIPointer(t, i, r, s, o)
	}

	function v() {
		y(), c !== l && (c = l, h(c.object))
	}

	function y() {
		l.geometry = null, l.program = null, l.wireframe = !1
	}
	return {
		setup: function(r, l, u, v, y) {
			let E = !1;
			if (s) {
				const t = function(t, i, r) {
					const s = !0 === r.wireframe;
					let l = o[t.id];
					void 0 === l && (l = {}, o[t.id] = l);
					let c = l[i.id];
					void 0 === c && (c = {}, l[i.id] = c);
					let h = c[s];
					void 0 === h && (h = d(n.isWebGL2 ? e.createVertexArray() : a.createVertexArrayOES()), c[s] = h);
					return h
				}(v, u, l);
				c !== t && (c = t, h(c.object)), E = function(e, t) {
					const i = c.attributes,
						n = e.attributes;
					let r = 0;
					for (const e in n) {
						const t = i[e],
							a = n[e];
						if (void 0 === t) return !0;
						if (t.attribute !== a) return !0;
						if (t.data !== a.data) return !0;
						r++
					}
					return c.attributesNum !== r || c.index !== t
				}(v, y), E && function(e, t) {
					const i = {},
						n = e.attributes;
					let r = 0;
					for (const e in n) {
						const t = n[e],
							a = {};
						a.attribute = t, t.data && (a.data = t.data), i[e] = a, r++
					}
					c.attributes = i, c.attributesNum = r, c.index = t
				}(v, y)
			} else {
				const e = !0 === l.wireframe;
				c.geometry === v.id && c.program === u.id && c.wireframe === e || (c.geometry = v.id, c.program = u.id, c.wireframe = e, E = !0)
			}!0 === r.isInstancedMesh && (E = !0), null !== y && i.update(y, 34963), E && (! function(r, a, s, o) {
				if (!1 === n.isWebGL2 && (r.isInstancedMesh || o.isInstancedBufferGeometry) && null === t.get("ANGLE_instanced_arrays")) return;
				p();
				const l = o.attributes,
					c = s.getAttributes(),
					h = a.defaultAttributeValues;
				for (const t in c) {
					const n = c[t];
					if (n >= 0) {
						const a = l[t];
						if (void 0 !== a) {
							const t = a.normalized,
								r = a.itemSize,
								s = i.get(a);
							if (void 0 === s) continue;
							const l = s.buffer,
								c = s.type,
								h = s.bytesPerElement;
							if (a.isInterleavedBufferAttribute) {
								const i = a.data,
									s = i.stride,
									u = a.offset;
								i && i.isInstancedInterleavedBuffer ? (A(n, i.meshPerAttribute), void 0 === o._maxInstanceCount && (o._maxInstanceCount = i.meshPerAttribute * i.count)) : m(n), e.bindBuffer(34962, l), f(n, r, c, t, s * h, u * h)
							} else a.isInstancedBufferAttribute ? (A(n, a.meshPerAttribute), void 0 === o._maxInstanceCount && (o._maxInstanceCount = a.meshPerAttribute * a.count)) : m(n), e.bindBuffer(34962, l), f(n, r, c, t, 0, 0)
						} else if ("instanceMatrix" === t) {
							const t = i.get(r.instanceMatrix);
							if (void 0 === t) continue;
							const a = t.buffer,
								s = t.type;
							A(n + 0, 1), A(n + 1, 1), A(n + 2, 1), A(n + 3, 1), e.bindBuffer(34962, a), e.vertexAttribPointer(n + 0, 4, s, !1, 64, 0), e.vertexAttribPointer(n + 1, 4, s, !1, 64, 16), e.vertexAttribPointer(n + 2, 4, s, !1, 64, 32), e.vertexAttribPointer(n + 3, 4, s, !1, 64, 48)
						} else if ("instanceColor" === t) {
							const t = i.get(r.instanceColor);
							if (void 0 === t) continue;
							const a = t.buffer,
								s = t.type;
							A(n, 1), e.bindBuffer(34962, a), e.vertexAttribPointer(n, 3, s, !1, 12, 0)
						} else if (void 0 !== h) {
							const i = h[t];
							if (void 0 !== i) switch (i.length) {
								case 2:
									e.vertexAttrib2fv(n, i);
									break;
								case 3:
									e.vertexAttrib3fv(n, i);
									break;
								case 4:
									e.vertexAttrib4fv(n, i);
									break;
								default:
									e.vertexAttrib1fv(n, i)
							}
						}
					}
				}
				g()
			}(r, l, u, v), null !== y && e.bindBuffer(34963, i.get(y).buffer))
		},
		reset: v,
		resetDefaultState: y,
		dispose: function() {
			v();
			for (const e in o) {
				const t = o[e];
				for (const e in t) {
					const i = t[e];
					for (const e in i) u(i[e].object), delete i[e];
					delete t[e]
				}
				delete o[e]
			}
		},
		releaseStatesOfGeometry: function(e) {
			if (void 0 === o[e.id]) return;
			const t = o[e.id];
			for (const e in t) {
				const i = t[e];
				for (const e in i) u(i[e].object), delete i[e];
				delete t[e]
			}
			delete o[e.id]
		},
		releaseStatesOfProgram: function(e) {
			for (const t in o) {
				const i = o[t];
				if (void 0 === i[e.id]) continue;
				const n = i[e.id];
				for (const e in n) u(n[e].object), delete n[e];
				delete i[e.id]
			}
		},
		initAttributes: p,
		enableAttribute: m,
		disableUnusedAttributes: g
	}
}

function WebGLBufferRenderer(e, t, i, n) {
	const r = n.isWebGL2;
	let a;
	this.setMode = function(e) {
		a = e
	}, this.render = function(t, n) {
		e.drawArrays(a, t, n), i.update(n, a, 1)
	}, this.renderInstances = function(n, s, o) {
		if (0 === o) return;
		let l, c;
		if (r) l = e, c = "drawArraysInstanced";
		else if (l = t.get("ANGLE_instanced_arrays"), c = "drawArraysInstancedANGLE", null === l) return void console.error("THREE.WebGLBufferRenderer: using THREE.InstancedBufferGeometry but hardware does not support extension ANGLE_instanced_arrays.");
		l[c](a, n, s, o), i.update(s, a, o)
	}
}

function WebGLCapabilities(e, t, i) {
	let n;

	function r(t) {
		if ("highp" === t) {
			if (e.getShaderPrecisionFormat(35633, 36338).precision > 0 && e.getShaderPrecisionFormat(35632, 36338).precision > 0) return "highp";
			t = "mediump"
		}
		return "mediump" === t && e.getShaderPrecisionFormat(35633, 36337).precision > 0 && e.getShaderPrecisionFormat(35632, 36337).precision > 0 ? "mediump" : "lowp"
	}
	const a = "undefined" != typeof WebGL2RenderingContext && e instanceof WebGL2RenderingContext || "undefined" != typeof WebGL2ComputeRenderingContext && e instanceof WebGL2ComputeRenderingContext;
	let s = void 0 !== i.precision ? i.precision : "highp";
	const o = r(s);
	o !== s && (console.warn("THREE.WebGLRenderer:", s, "not supported, using", o, "instead."), s = o);
	const l = !0 === i.logarithmicDepthBuffer,
		c = e.getParameter(34930),
		h = e.getParameter(35660),
		u = e.getParameter(3379),
		d = e.getParameter(34076),
		p = e.getParameter(34921),
		m = e.getParameter(36347),
		A = e.getParameter(36348),
		g = e.getParameter(36349),
		f = h > 0,
		v = a || t.has("OES_texture_float");
	return {
		isWebGL2: a,
		getMaxAnisotropy: function() {
			if (void 0 !== n) return n;
			if (!0 === t.has("EXT_texture_filter_anisotropic")) {
				const i = t.get("EXT_texture_filter_anisotropic");
				n = e.getParameter(i.MAX_TEXTURE_MAX_ANISOTROPY_EXT)
			} else n = 0;
			return n
		},
		getMaxPrecision: r,
		precision: s,
		logarithmicDepthBuffer: l,
		maxTextures: c,
		maxVertexTextures: h,
		maxTextureSize: u,
		maxCubemapSize: d,
		maxAttributes: p,
		maxVertexUniforms: m,
		maxVaryings: A,
		maxFragmentUniforms: g,
		vertexTextures: f,
		floatFragmentTextures: v,
		floatVertexTextures: f && v,
		maxSamples: a ? e.getParameter(36183) : 0
	}
}

function WebGLClipping(e) {
	const t = this;
	let i = null,
		n = 0,
		r = !1,
		a = !1;
	const s = new Plane,
		o = new Matrix3,
		l = {
			value: null,
			needsUpdate: !1
		};

	function c() {
		l.value !== i && (l.value = i, l.needsUpdate = n > 0), t.numPlanes = n, t.numIntersection = 0
	}

	function h(e, i, n, r) {
		const a = null !== e ? e.length : 0;
		let c = null;
		if (0 !== a) {
			if (c = l.value, !0 !== r || null === c) {
				const t = n + 4 * a,
					r = i.matrixWorldInverse;
				o.getNormalMatrix(r), (null === c || c.length < t) && (c = new Float32Array(t));
				for (let t = 0, i = n; t !== a; ++t, i += 4) s.copy(e[t]).applyMatrix4(r, o), s.normal.toArray(c, i), c[i + 3] = s.constant
			}
			l.value = c, l.needsUpdate = !0
		}
		return t.numPlanes = a, t.numIntersection = 0, c
	}
	this.uniform = l, this.numPlanes = 0, this.numIntersection = 0, this.init = function(e, t, a) {
		const s = 0 !== e.length || t || 0 !== n || r;
		return r = t, i = h(e, a, 0), n = e.length, s
	}, this.beginShadows = function() {
		a = !0, h(null)
	}, this.endShadows = function() {
		a = !1, c()
	}, this.setState = function(t, s, o) {
		const u = t.clippingPlanes,
			d = t.clipIntersection,
			p = t.clipShadows,
			m = e.get(t);
		if (!r || null === u || 0 === u.length || a && !p) a ? h(null) : c();
		else {
			const e = a ? 0 : n,
				t = 4 * e;
			let r = m.clippingState || null;
			l.value = r, r = h(u, s, t, o);
			for (let e = 0; e !== t; ++e) r[e] = i[e];
			m.clippingState = r, this.numIntersection = d ? this.numPlanes : 0, this.numPlanes += e
		}
	}
}

function WebGLCubeMaps(e) {
	let t = new WeakMap;

	function i(e, t) {
		return 303 === t ? e.mapping = 301 : 304 === t && (e.mapping = 302), e
	}

	function n(e) {
		const i = e.target;
		i.removeEventListener("dispose", n);
		const r = t.get(i);
		void 0 !== r && (t.delete(i), r.dispose())
	}
	return {
		get: function(r) {
			if (r && r.isTexture) {
				const a = r.mapping;
				if (303 === a || 304 === a) {
					if (t.has(r)) {
						return i(t.get(r).texture, r.mapping)
					} {
						const a = r.image;
						if (a && a.height > 0) {
							const s = e.getRenderTarget(),
								o = new WebGLCubeRenderTarget(a.height / 2);
							return o.fromEquirectangularTexture(e, r), t.set(r, o), e.setRenderTarget(s), r.addEventListener("dispose", n), i(o.texture, r.mapping)
						}
						return null
					}
				}
			}
			return r
		},
		dispose: function() {
			t = new WeakMap
		}
	}
}

function WebGLExtensions(e) {
	const t = {};

	function i(i) {
		if (void 0 !== t[i]) return t[i];
		let n;
		switch (i) {
			case "WEBGL_depth_texture":
				n = e.getExtension("WEBGL_depth_texture") || e.getExtension("MOZ_WEBGL_depth_texture") || e.getExtension("WEBKIT_WEBGL_depth_texture");
				break;
			case "EXT_texture_filter_anisotropic":
				n = e.getExtension("EXT_texture_filter_anisotropic") || e.getExtension("MOZ_EXT_texture_filter_anisotropic") || e.getExtension("WEBKIT_EXT_texture_filter_anisotropic");
				break;
			case "WEBGL_compressed_texture_s3tc":
				n = e.getExtension("WEBGL_compressed_texture_s3tc") || e.getExtension("MOZ_WEBGL_compressed_texture_s3tc") || e.getExtension("WEBKIT_WEBGL_compressed_texture_s3tc");
				break;
			case "WEBGL_compressed_texture_pvrtc":
				n = e.getExtension("WEBGL_compressed_texture_pvrtc") || e.getExtension("WEBKIT_WEBGL_compressed_texture_pvrtc");
				break;
			default:
				n = e.getExtension(i)
		}
		return t[i] = n, n
	}
	return {
		has: function(e) {
			return null !== i(e)
		},
		init: function(e) {
			e.isWebGL2 ? i("EXT_color_buffer_float") : (i("WEBGL_depth_texture"), i("OES_texture_float"), i("OES_texture_half_float"), i("OES_texture_half_float_linear"), i("OES_standard_derivatives"), i("OES_element_index_uint"), i("OES_vertex_array_object"), i("ANGLE_instanced_arrays")), i("OES_texture_float_linear"), i("EXT_color_buffer_half_float")
		},
		get: function(e) {
			const t = i(e);
			return null === t && console.warn("THREE.WebGLRenderer: " + e + " extension not supported."), t
		}
	}
}

function WebGLGeometries(e, t, i, n) {
	const r = {},
		a = new WeakMap;

	function s(e) {
		const o = e.target;
		null !== o.index && t.remove(o.index);
		for (const e in o.attributes) t.remove(o.attributes[e]);
		o.removeEventListener("dispose", s), delete r[o.id];
		const l = a.get(o);
		l && (t.remove(l), a.delete(o)), n.releaseStatesOfGeometry(o), !0 === o.isInstancedBufferGeometry && delete o._maxInstanceCount, i.memory.geometries--
	}

	function o(e) {
		const i = [],
			n = e.index,
			r = e.attributes.position;
		let s = 0;
		if (null !== n) {
			const e = n.array;
			s = n.version;
			for (let t = 0, n = e.length; t < n; t += 3) {
				const n = e[t + 0],
					r = e[t + 1],
					a = e[t + 2];
				i.push(n, r, r, a, a, n)
			}
		} else {
			const e = r.array;
			s = r.version;
			for (let t = 0, n = e.length / 3 - 1; t < n; t += 3) {
				const e = t + 0,
					n = t + 1,
					r = t + 2;
				i.push(e, n, n, r, r, e)
			}
		}
		const o = new(arrayMax(i) > 65535 ? Uint32BufferAttribute : Uint16BufferAttribute)(i, 1);
		o.version = s;
		const l = a.get(e);
		l && t.remove(l), a.set(e, o)
	}
	return {
		get: function(e, t) {
			return !0 === r[t.id] || (t.addEventListener("dispose", s), r[t.id] = !0, i.memory.geometries++), t
		},
		update: function(e) {
			const i = e.attributes;
			for (const e in i) t.update(i[e], 34962);
			const n = e.morphAttributes;
			for (const e in n) {
				const i = n[e];
				for (let e = 0, n = i.length; e < n; e++) t.update(i[e], 34962)
			}
		},
		getWireframeAttribute: function(e) {
			const t = a.get(e);
			if (t) {
				const i = e.index;
				null !== i && t.version < i.version && o(e)
			} else o(e);
			return a.get(e)
		}
	}
}

function WebGLIndexedBufferRenderer(e, t, i, n) {
	const r = n.isWebGL2;
	let a, s, o;
	this.setMode = function(e) {
		a = e
	}, this.setIndex = function(e) {
		s = e.type, o = e.bytesPerElement
	}, this.render = function(t, n) {
		e.drawElements(a, n, s, t * o), i.update(n, a, 1)
	}, this.renderInstances = function(n, l, c) {
		if (0 === c) return;
		let h, u;
		if (r) h = e, u = "drawElementsInstanced";
		else if (h = t.get("ANGLE_instanced_arrays"), u = "drawElementsInstancedANGLE", null === h) return void console.error("THREE.WebGLIndexedBufferRenderer: using THREE.InstancedBufferGeometry but hardware does not support extension ANGLE_instanced_arrays.");
		h[u](a, l, s, n * o, c), i.update(l, a, c)
	}
}

function WebGLInfo(e) {
	const t = {
		frame: 0,
		calls: 0,
		triangles: 0,
		points: 0,
		lines: 0
	};
	return {
		memory: {
			geometries: 0,
			textures: 0
		},
		render: t,
		programs: null,
		autoReset: !0,
		reset: function() {
			t.frame++, t.calls = 0, t.triangles = 0, t.points = 0, t.lines = 0
		},
		update: function(e, i, n) {
			switch (t.calls++, i) {
				case 4:
					t.triangles += n * (e / 3);
					break;
				case 1:
					t.lines += n * (e / 2);
					break;
				case 3:
					t.lines += n * (e - 1);
					break;
				case 2:
					t.lines += n * e;
					break;
				case 0:
					t.points += n * e;
					break;
				default:
					console.error("THREE.WebGLInfo: Unknown draw mode:", i)
			}
		}
	}
}

function numericalSort(e, t) {
	return e[0] - t[0]
}

function absNumericalSort(e, t) {
	return Math.abs(t[1]) - Math.abs(e[1])
}

function WebGLMorphtargets(e) {
	const t = {},
		i = new Float32Array(8),
		n = [];
	for (let e = 0; e < 8; e++) n[e] = [e, 0];
	return {
		update: function(r, a, s, o) {
			const l = r.morphTargetInfluences,
				c = void 0 === l ? 0 : l.length;
			let h = t[a.id];
			if (void 0 === h) {
				h = [];
				for (let e = 0; e < c; e++) h[e] = [e, 0];
				t[a.id] = h
			}
			for (let e = 0; e < c; e++) {
				const t = h[e];
				t[0] = e, t[1] = l[e]
			}
			h.sort(absNumericalSort);
			for (let e = 0; e < 8; e++) e < c && h[e][1] ? (n[e][0] = h[e][0], n[e][1] = h[e][1]) : (n[e][0] = Number.MAX_SAFE_INTEGER, n[e][1] = 0);
			n.sort(numericalSort);
			const u = s.morphTargets && a.morphAttributes.position,
				d = s.morphNormals && a.morphAttributes.normal;
			let p = 0;
			for (let e = 0; e < 8; e++) {
				const t = n[e],
					r = t[0],
					s = t[1];
				r !== Number.MAX_SAFE_INTEGER && s ? (u && a.getAttribute("morphTarget" + e) !== u[r] && a.setAttribute("morphTarget" + e, u[r]), d && a.getAttribute("morphNormal" + e) !== d[r] && a.setAttribute("morphNormal" + e, d[r]), i[e] = s, p += s) : (u && !0 === a.hasAttribute("morphTarget" + e) && a.deleteAttribute("morphTarget" + e), d && !0 === a.hasAttribute("morphNormal" + e) && a.deleteAttribute("morphNormal" + e), i[e] = 0)
			}
			const m = a.morphTargetsRelative ? 1 : 1 - p;
			o.getUniforms().setValue(e, "morphTargetBaseInfluence", m), o.getUniforms().setValue(e, "morphTargetInfluences", i)
		}
	}
}

function WebGLObjects(e, t, i, n) {
	let r = new WeakMap;

	function a(e) {
		const t = e.target;
		t.removeEventListener("dispose", a), i.remove(t.instanceMatrix), null !== t.instanceColor && i.remove(t.instanceColor)
	}
	return {
		update: function(e) {
			const s = n.render.frame,
				o = e.geometry,
				l = t.get(e, o);
			return r.get(l) !== s && (t.update(l), r.set(l, s)), e.isInstancedMesh && (!1 === e.hasEventListener("dispose", a) && e.addEventListener("dispose", a), i.update(e.instanceMatrix, 34962), null !== e.instanceColor && i.update(e.instanceColor, 34962)), l
		},
		dispose: function() {
			r = new WeakMap
		}
	}
}
ShaderLib.physical = {
	uniforms: mergeUniforms([ShaderLib.standard.uniforms, {
		clearcoat: {
			value: 0
		},
		clearcoatMap: {
			value: null
		},
		clearcoatRoughness: {
			value: 0
		},
		clearcoatRoughnessMap: {
			value: null
		},
		clearcoatNormalScale: {
			value: new Vector2(1, 1)
		},
		clearcoatNormalMap: {
			value: null
		},
		sheen: {
			value: new Color(0)
		},
		transmission: {
			value: 0
		},
		transmissionMap: {
			value: null
		}
	}]),
	vertexShader: ShaderChunk.meshphysical_vert,
	fragmentShader: ShaderChunk.meshphysical_frag
};
class DataTexture2DArray extends Texture$1 {
	constructor(e = null, t = 1, i = 1, n = 1) {
		super(null), this.image = {
			data: e,
			width: t,
			height: i,
			depth: n
		}, this.magFilter = 1003, this.minFilter = 1003, this.wrapR = 1001, this.generateMipmaps = !1, this.flipY = !1, this.needsUpdate = !0
	}
}
DataTexture2DArray.prototype.isDataTexture2DArray = !0;
class DataTexture3D extends Texture$1 {
	constructor(e = null, t = 1, i = 1, n = 1) {
		super(null), this.image = {
			data: e,
			width: t,
			height: i,
			depth: n
		}, this.magFilter = 1003, this.minFilter = 1003, this.wrapR = 1001, this.generateMipmaps = !1, this.flipY = !1, this.needsUpdate = !0
	}
}
DataTexture3D.prototype.isDataTexture3D = !0;
const emptyTexture = new Texture$1,
	emptyTexture2dArray = new DataTexture2DArray,
	emptyTexture3d = new DataTexture3D,
	emptyCubeTexture = new CubeTexture,
	arrayCacheF32 = [],
	arrayCacheI32 = [],
	mat4array = new Float32Array(16),
	mat3array = new Float32Array(9),
	mat2array = new Float32Array(4);

function flatten(e, t, i) {
	const n = e[0];
	if (n <= 0 || n > 0) return e;
	const r = t * i;
	let a = arrayCacheF32[r];
	if (void 0 === a && (a = new Float32Array(r), arrayCacheF32[r] = a), 0 !== t) {
		n.toArray(a, 0);
		for (let n = 1, r = 0; n !== t; ++n) r += i, e[n].toArray(a, r)
	}
	return a
}

function arraysEqual(e, t) {
	if (e.length !== t.length) return !1;
	for (let i = 0, n = e.length; i < n; i++)
		if (e[i] !== t[i]) return !1;
	return !0
}

function copyArray(e, t) {
	for (let i = 0, n = t.length; i < n; i++) e[i] = t[i]
}

function allocTexUnits(e, t) {
	let i = arrayCacheI32[t];
	void 0 === i && (i = new Int32Array(t), arrayCacheI32[t] = i);
	for (let n = 0; n !== t; ++n) i[n] = e.allocateTextureUnit();
	return i
}

function setValueV1f(e, t) {
	const i = this.cache;
	i[0] !== t && (e.uniform1f(this.addr, t), i[0] = t)
}

function setValueV2f(e, t) {
	const i = this.cache;
	if (void 0 !== t.x) i[0] === t.x && i[1] === t.y || (e.uniform2f(this.addr, t.x, t.y), i[0] = t.x, i[1] = t.y);
	else {
		if (arraysEqual(i, t)) return;
		e.uniform2fv(this.addr, t), copyArray(i, t)
	}
}

function setValueV3f(e, t) {
	const i = this.cache;
	if (void 0 !== t.x) i[0] === t.x && i[1] === t.y && i[2] === t.z || (e.uniform3f(this.addr, t.x, t.y, t.z), i[0] = t.x, i[1] = t.y, i[2] = t.z);
	else if (void 0 !== t.r) i[0] === t.r && i[1] === t.g && i[2] === t.b || (e.uniform3f(this.addr, t.r, t.g, t.b), i[0] = t.r, i[1] = t.g, i[2] = t.b);
	else {
		if (arraysEqual(i, t)) return;
		e.uniform3fv(this.addr, t), copyArray(i, t)
	}
}

function setValueV4f(e, t) {
	const i = this.cache;
	if (void 0 !== t.x) i[0] === t.x && i[1] === t.y && i[2] === t.z && i[3] === t.w || (e.uniform4f(this.addr, t.x, t.y, t.z, t.w), i[0] = t.x, i[1] = t.y, i[2] = t.z, i[3] = t.w);
	else {
		if (arraysEqual(i, t)) return;
		e.uniform4fv(this.addr, t), copyArray(i, t)
	}
}

function setValueM2(e, t) {
	const i = this.cache,
		n = t.elements;
	if (void 0 === n) {
		if (arraysEqual(i, t)) return;
		e.uniformMatrix2fv(this.addr, !1, t), copyArray(i, t)
	} else {
		if (arraysEqual(i, n)) return;
		mat2array.set(n), e.uniformMatrix2fv(this.addr, !1, mat2array), copyArray(i, n)
	}
}

function setValueM3(e, t) {
	const i = this.cache,
		n = t.elements;
	if (void 0 === n) {
		if (arraysEqual(i, t)) return;
		e.uniformMatrix3fv(this.addr, !1, t), copyArray(i, t)
	} else {
		if (arraysEqual(i, n)) return;
		mat3array.set(n), e.uniformMatrix3fv(this.addr, !1, mat3array), copyArray(i, n)
	}
}

function setValueM4(e, t) {
	const i = this.cache,
		n = t.elements;
	if (void 0 === n) {
		if (arraysEqual(i, t)) return;
		e.uniformMatrix4fv(this.addr, !1, t), copyArray(i, t)
	} else {
		if (arraysEqual(i, n)) return;
		mat4array.set(n), e.uniformMatrix4fv(this.addr, !1, mat4array), copyArray(i, n)
	}
}

function setValueT1(e, t, i) {
	const n = this.cache,
		r = i.allocateTextureUnit();
	n[0] !== r && (e.uniform1i(this.addr, r), n[0] = r), i.safeSetTexture2D(t || emptyTexture, r)
}

function setValueT2DArray1(e, t, i) {
	const n = this.cache,
		r = i.allocateTextureUnit();
	n[0] !== r && (e.uniform1i(this.addr, r), n[0] = r), i.setTexture2DArray(t || emptyTexture2dArray, r)
}

function setValueT3D1(e, t, i) {
	const n = this.cache,
		r = i.allocateTextureUnit();
	n[0] !== r && (e.uniform1i(this.addr, r), n[0] = r), i.setTexture3D(t || emptyTexture3d, r)
}

function setValueT6(e, t, i) {
	const n = this.cache,
		r = i.allocateTextureUnit();
	n[0] !== r && (e.uniform1i(this.addr, r), n[0] = r), i.safeSetTextureCube(t || emptyCubeTexture, r)
}

function setValueV1i(e, t) {
	const i = this.cache;
	i[0] !== t && (e.uniform1i(this.addr, t), i[0] = t)
}

function setValueV2i(e, t) {
	const i = this.cache;
	arraysEqual(i, t) || (e.uniform2iv(this.addr, t), copyArray(i, t))
}

function setValueV3i(e, t) {
	const i = this.cache;
	arraysEqual(i, t) || (e.uniform3iv(this.addr, t), copyArray(i, t))
}

function setValueV4i(e, t) {
	const i = this.cache;
	arraysEqual(i, t) || (e.uniform4iv(this.addr, t), copyArray(i, t))
}

function setValueV1ui(e, t) {
	const i = this.cache;
	i[0] !== t && (e.uniform1ui(this.addr, t), i[0] = t)
}

function getSingularSetter(e) {
	switch (e) {
		case 5126:
			return setValueV1f;
		case 35664:
			return setValueV2f;
		case 35665:
			return setValueV3f;
		case 35666:
			return setValueV4f;
		case 35674:
			return setValueM2;
		case 35675:
			return setValueM3;
		case 35676:
			return setValueM4;
		case 5124:
		case 35670:
			return setValueV1i;
		case 35667:
		case 35671:
			return setValueV2i;
		case 35668:
		case 35672:
			return setValueV3i;
		case 35669:
		case 35673:
			return setValueV4i;
		case 5125:
			return setValueV1ui;
		case 35678:
		case 36198:
		case 36298:
		case 36306:
		case 35682:
			return setValueT1;
		case 35679:
		case 36299:
		case 36307:
			return setValueT3D1;
		case 35680:
		case 36300:
		case 36308:
		case 36293:
			return setValueT6;
		case 36289:
		case 36303:
		case 36311:
		case 36292:
			return setValueT2DArray1
	}
}

function setValueV1fArray(e, t) {
	e.uniform1fv(this.addr, t)
}

function setValueV1iArray(e, t) {
	e.uniform1iv(this.addr, t)
}

function setValueV2iArray(e, t) {
	e.uniform2iv(this.addr, t)
}

function setValueV3iArray(e, t) {
	e.uniform3iv(this.addr, t)
}

function setValueV4iArray(e, t) {
	e.uniform4iv(this.addr, t)
}

function setValueV2fArray(e, t) {
	const i = flatten(t, this.size, 2);
	e.uniform2fv(this.addr, i)
}

function setValueV3fArray(e, t) {
	const i = flatten(t, this.size, 3);
	e.uniform3fv(this.addr, i)
}

function setValueV4fArray(e, t) {
	const i = flatten(t, this.size, 4);
	e.uniform4fv(this.addr, i)
}

function setValueM2Array(e, t) {
	const i = flatten(t, this.size, 4);
	e.uniformMatrix2fv(this.addr, !1, i)
}

function setValueM3Array(e, t) {
	const i = flatten(t, this.size, 9);
	e.uniformMatrix3fv(this.addr, !1, i)
}

function setValueM4Array(e, t) {
	const i = flatten(t, this.size, 16);
	e.uniformMatrix4fv(this.addr, !1, i)
}

function setValueT1Array(e, t, i) {
	const n = t.length,
		r = allocTexUnits(i, n);
	e.uniform1iv(this.addr, r);
	for (let e = 0; e !== n; ++e) i.safeSetTexture2D(t[e] || emptyTexture, r[e])
}

function setValueT6Array(e, t, i) {
	const n = t.length,
		r = allocTexUnits(i, n);
	e.uniform1iv(this.addr, r);
	for (let e = 0; e !== n; ++e) i.safeSetTextureCube(t[e] || emptyCubeTexture, r[e])
}

function getPureArraySetter(e) {
	switch (e) {
		case 5126:
			return setValueV1fArray;
		case 35664:
			return setValueV2fArray;
		case 35665:
			return setValueV3fArray;
		case 35666:
			return setValueV4fArray;
		case 35674:
			return setValueM2Array;
		case 35675:
			return setValueM3Array;
		case 35676:
			return setValueM4Array;
		case 5124:
		case 35670:
			return setValueV1iArray;
		case 35667:
		case 35671:
			return setValueV2iArray;
		case 35668:
		case 35672:
			return setValueV3iArray;
		case 35669:
		case 35673:
			return setValueV4iArray;
		case 35678:
		case 36198:
		case 36298:
		case 36306:
		case 35682:
			return setValueT1Array;
		case 35680:
		case 36300:
		case 36308:
		case 36293:
			return setValueT6Array
	}
}

function SingleUniform(e, t, i) {
	this.id = e, this.addr = i, this.cache = [], this.setValue = getSingularSetter(t.type)
}

function PureArrayUniform(e, t, i) {
	this.id = e, this.addr = i, this.cache = [], this.size = t.size, this.setValue = getPureArraySetter(t.type)
}

function StructuredUniform(e) {
	this.id = e, this.seq = [], this.map = {}
}
PureArrayUniform.prototype.updateCache = function(e) {
	const t = this.cache;
	e instanceof Float32Array && t.length !== e.length && (this.cache = new Float32Array(e.length)), copyArray(t, e)
}, StructuredUniform.prototype.setValue = function(e, t, i) {
	const n = this.seq;
	for (let r = 0, a = n.length; r !== a; ++r) {
		const a = n[r];
		a.setValue(e, t[a.id], i)
	}
};
const RePathPart = /(\w+)(\])?(\[|\.)?/g;

function addUniform(e, t) {
	e.seq.push(t), e.map[t.id] = t
}

function parseUniform(e, t, i) {
	const n = e.name,
		r = n.length;
	for (RePathPart.lastIndex = 0;;) {
		const a = RePathPart.exec(n),
			s = RePathPart.lastIndex;
		let o = a[1];
		const l = "]" === a[2],
			c = a[3];
		if (l && (o |= 0), void 0 === c || "[" === c && s + 2 === r) {
			addUniform(i, void 0 === c ? new SingleUniform(o, e, t) : new PureArrayUniform(o, e, t));
			break
		} {
			let e = i.map[o];
			void 0 === e && (e = new StructuredUniform(o), addUniform(i, e)), i = e
		}
	}
}

function WebGLUniforms(e, t) {
	this.seq = [], this.map = {};
	const i = e.getProgramParameter(t, 35718);
	for (let n = 0; n < i; ++n) {
		const i = e.getActiveUniform(t, n);
		parseUniform(i, e.getUniformLocation(t, i.name), this)
	}
}

function WebGLShader(e, t, i) {
	const n = e.createShader(t);
	return e.shaderSource(n, i), e.compileShader(n), n
}
WebGLUniforms.prototype.setValue = function(e, t, i, n) {
	const r = this.map[t];
	void 0 !== r && r.setValue(e, i, n)
}, WebGLUniforms.prototype.setOptional = function(e, t, i) {
	const n = t[i];
	void 0 !== n && this.setValue(e, i, n)
}, WebGLUniforms.upload = function(e, t, i, n) {
	for (let r = 0, a = t.length; r !== a; ++r) {
		const a = t[r],
			s = i[a.id];
		!1 !== s.needsUpdate && a.setValue(e, s.value, n)
	}
}, WebGLUniforms.seqWithValue = function(e, t) {
	const i = [];
	for (let n = 0, r = e.length; n !== r; ++n) {
		const r = e[n];
		r.id in t && i.push(r)
	}
	return i
};
let programIdCount = 0;

function addLineNumbers(e) {
	const t = e.split("\n");
	for (let e = 0; e < t.length; e++) t[e] = e + 1 + ": " + t[e];
	return t.join("\n")
}

function getEncodingComponents(e) {
	switch (e) {
		case 3e3:
			return ["Linear", "( value )"];
		case 3001:
			return ["sRGB", "( value )"];
		case 3002:
			return ["RGBE", "( value )"];
		case 3004:
			return ["RGBM", "( value, 7.0 )"];
		case 3005:
			return ["RGBM", "( value, 16.0 )"];
		case 3006:
			return ["RGBD", "( value, 256.0 )"];
		case 3007:
			return ["Gamma", "( value, float( GAMMA_FACTOR ) )"];
		case 3003:
			return ["LogLuv", "( value )"];
		default:
			return console.warn("THREE.WebGLProgram: Unsupported encoding:", e), ["Linear", "( value )"]
	}
}

function getShaderErrors(e, t, i) {
	const n = e.getShaderParameter(t, 35713),
		r = e.getShaderInfoLog(t).trim();
	if (n && "" === r) return "";
	return "THREE.WebGLShader: gl.getShaderInfoLog() " + i + "\n" + r + addLineNumbers(e.getShaderSource(t))
}

function getTexelDecodingFunction(e, t) {
	const i = getEncodingComponents(t);
	return "vec4 " + e + "( vec4 value ) { return " + i[0] + "ToLinear" + i[1] + "; }"
}

function getTexelEncodingFunction(e, t) {
	const i = getEncodingComponents(t);
	return "vec4 " + e + "( vec4 value ) { return LinearTo" + i[0] + i[1] + "; }"
}

function getToneMappingFunction(e, t) {
	let i;
	switch (t) {
		case 1:
			i = "Linear";
			break;
		case 2:
			i = "Reinhard";
			break;
		case 3:
			i = "OptimizedCineon";
			break;
		case 4:
			i = "ACESFilmic";
			break;
		case 5:
			i = "Custom";
			break;
		default:
			console.warn("THREE.WebGLProgram: Unsupported toneMapping:", t), i = "Linear"
	}
	return "vec3 " + e + "( vec3 color ) { return " + i + "ToneMapping( color ); }"
}

function generateExtensions(e) {
	return [e.extensionDerivatives || e.envMapCubeUV || e.bumpMap || e.tangentSpaceNormalMap || e.clearcoatNormalMap || e.flatShading || "physical" === e.shaderID ? "#extension GL_OES_standard_derivatives : enable" : "", (e.extensionFragDepth || e.logarithmicDepthBuffer) && e.rendererExtensionFragDepth ? "#extension GL_EXT_frag_depth : enable" : "", e.extensionDrawBuffers && e.rendererExtensionDrawBuffers ? "#extension GL_EXT_draw_buffers : require" : "", (e.extensionShaderTextureLOD || e.envMap) && e.rendererExtensionShaderTextureLod ? "#extension GL_EXT_shader_texture_lod : enable" : ""].filter(filterEmptyLine).join("\n")
}

function generateDefines(e) {
	const t = [];
	for (const i in e) {
		const n = e[i];
		!1 !== n && t.push("#define " + i + " " + n)
	}
	return t.join("\n")
}

function fetchAttributeLocations(e, t) {
	const i = {},
		n = e.getProgramParameter(t, 35721);
	for (let r = 0; r < n; r++) {
		const n = e.getActiveAttrib(t, r).name;
		i[n] = e.getAttribLocation(t, n)
	}
	return i
}

function filterEmptyLine(e) {
	return "" !== e
}

function replaceLightNums(e, t) {
	return e.replace(/NUM_DIR_LIGHTS/g, t.numDirLights).replace(/NUM_SPOT_LIGHTS/g, t.numSpotLights).replace(/NUM_RECT_AREA_LIGHTS/g, t.numRectAreaLights).replace(/NUM_POINT_LIGHTS/g, t.numPointLights).replace(/NUM_HEMI_LIGHTS/g, t.numHemiLights).replace(/NUM_DIR_LIGHT_SHADOWS/g, t.numDirLightShadows).replace(/NUM_SPOT_LIGHT_SHADOWS/g, t.numSpotLightShadows).replace(/NUM_POINT_LIGHT_SHADOWS/g, t.numPointLightShadows)
}

function replaceClippingPlaneNums(e, t) {
	return e.replace(/NUM_CLIPPING_PLANES/g, t.numClippingPlanes).replace(/UNION_CLIPPING_PLANES/g, t.numClippingPlanes - t.numClipIntersection)
}
const includePattern = /^[ \t]*#include +<([\w\d./]+)>/gm;

function resolveIncludes(e) {
	return e.replace(includePattern, includeReplacer)
}

function includeReplacer(e, t) {
	const i = ShaderChunk[t];
	if (void 0 === i) throw new Error("Can not resolve #include <" + t + ">");
	return resolveIncludes(i)
}
const deprecatedUnrollLoopPattern = /#pragma unroll_loop[\s]+?for \( int i \= (\d+)\; i < (\d+)\; i \+\+ \) \{([\s\S]+?)(?=\})\}/g,
	unrollLoopPattern = /#pragma unroll_loop_start\s+for\s*\(\s*int\s+i\s*=\s*(\d+)\s*;\s*i\s*<\s*(\d+)\s*;\s*i\s*\+\+\s*\)\s*{([\s\S]+?)}\s+#pragma unroll_loop_end/g;

function unrollLoops(e) {
	return e.replace(unrollLoopPattern, loopReplacer).replace(deprecatedUnrollLoopPattern, deprecatedLoopReplacer)
}

function deprecatedLoopReplacer(e, t, i, n) {
	return console.warn("WebGLProgram: #pragma unroll_loop shader syntax is deprecated. Please use #pragma unroll_loop_start syntax instead."), loopReplacer(e, t, i, n)
}

function loopReplacer(e, t, i, n) {
	let r = "";
	for (let e = parseInt(t); e < parseInt(i); e++) r += n.replace(/\[\s*i\s*\]/g, "[ " + e + " ]").replace(/UNROLLED_LOOP_INDEX/g, e);
	return r
}

function generatePrecision(e) {
	let t = "precision " + e.precision + " float;\nprecision " + e.precision + " int;";
	return "highp" === e.precision ? t += "\n#define HIGH_PRECISION" : "mediump" === e.precision ? t += "\n#define MEDIUM_PRECISION" : "lowp" === e.precision && (t += "\n#define LOW_PRECISION"), t
}

function generateShadowMapTypeDefine(e) {
	let t = "SHADOWMAP_TYPE_BASIC";
	return 1 === e.shadowMapType ? t = "SHADOWMAP_TYPE_PCF" : 2 === e.shadowMapType ? t = "SHADOWMAP_TYPE_PCF_SOFT" : 3 === e.shadowMapType && (t = "SHADOWMAP_TYPE_VSM"), t
}

function generateEnvMapTypeDefine(e) {
	let t = "ENVMAP_TYPE_CUBE";
	if (e.envMap) switch (e.envMapMode) {
		case 301:
		case 302:
			t = "ENVMAP_TYPE_CUBE";
			break;
		case 306:
		case 307:
			t = "ENVMAP_TYPE_CUBE_UV"
	}
	return t
}

function generateEnvMapModeDefine(e) {
	let t = "ENVMAP_MODE_REFLECTION";
	if (e.envMap) switch (e.envMapMode) {
		case 302:
		case 307:
			t = "ENVMAP_MODE_REFRACTION"
	}
	return t
}

function generateEnvMapBlendingDefine(e) {
	let t = "ENVMAP_BLENDING_NONE";
	if (e.envMap) switch (e.combine) {
		case 0:
			t = "ENVMAP_BLENDING_MULTIPLY";
			break;
		case 1:
			t = "ENVMAP_BLENDING_MIX";
			break;
		case 2:
			t = "ENVMAP_BLENDING_ADD"
	}
	return t
}

function WebGLProgram(e, t, i, n) {
	const r = e.getContext(),
		a = i.defines;
	let s = i.vertexShader,
		o = i.fragmentShader;
	const l = generateShadowMapTypeDefine(i),
		c = generateEnvMapTypeDefine(i),
		h = generateEnvMapModeDefine(i),
		u = generateEnvMapBlendingDefine(i),
		d = e.gammaFactor > 0 ? e.gammaFactor : 1,
		p = i.isWebGL2 ? "" : generateExtensions(i),
		m = generateDefines(a),
		A = r.createProgram();
	let g, f, v = i.glslVersion ? "#version " + i.glslVersion + "\n" : "";
	i.isRawShaderMaterial ? (g = [m].filter(filterEmptyLine).join("\n"), g.length > 0 && (g += "\n"), f = [p, m].filter(filterEmptyLine).join("\n"), f.length > 0 && (f += "\n")) : (g = [generatePrecision(i), "#define SHADER_NAME " + i.shaderName, m, i.instancing ? "#define USE_INSTANCING" : "", i.instancingColor ? "#define USE_INSTANCING_COLOR" : "", i.supportsVertexTextures ? "#define VERTEX_TEXTURES" : "", "#define GAMMA_FACTOR " + d, "#define MAX_BONES " + i.maxBones, i.useFog && i.fog ? "#define USE_FOG" : "", i.useFog && i.fogExp2 ? "#define FOG_EXP2" : "", i.map ? "#define USE_MAP" : "", i.envMap ? "#define USE_ENVMAP" : "", i.envMap ? "#define " + h : "", i.lightMap ? "#define USE_LIGHTMAP" : "", i.aoMap ? "#define USE_AOMAP" : "", i.emissiveMap ? "#define USE_EMISSIVEMAP" : "", i.bumpMap ? "#define USE_BUMPMAP" : "", i.normalMap ? "#define USE_NORMALMAP" : "", i.normalMap && i.objectSpaceNormalMap ? "#define OBJECTSPACE_NORMALMAP" : "", i.normalMap && i.tangentSpaceNormalMap ? "#define TANGENTSPACE_NORMALMAP" : "", i.clearcoatMap ? "#define USE_CLEARCOATMAP" : "", i.clearcoatRoughnessMap ? "#define USE_CLEARCOAT_ROUGHNESSMAP" : "", i.clearcoatNormalMap ? "#define USE_CLEARCOAT_NORMALMAP" : "", i.displacementMap && i.supportsVertexTextures ? "#define USE_DISPLACEMENTMAP" : "", i.specularMap ? "#define USE_SPECULARMAP" : "", i.roughnessMap ? "#define USE_ROUGHNESSMAP" : "", i.metalnessMap ? "#define USE_METALNESSMAP" : "", i.alphaMap ? "#define USE_ALPHAMAP" : "", i.transmissionMap ? "#define USE_TRANSMISSIONMAP" : "", i.vertexTangents ? "#define USE_TANGENT" : "", i.vertexColors ? "#define USE_COLOR" : "", i.vertexUvs ? "#define USE_UV" : "", i.uvsVertexOnly ? "#define UVS_VERTEX_ONLY" : "", i.flatShading ? "#define FLAT_SHADED" : "", i.skinning ? "#define USE_SKINNING" : "", i.useVertexTexture ? "#define BONE_TEXTURE" : "", i.morphTargets ? "#define USE_MORPHTARGETS" : "", i.morphNormals && !1 === i.flatShading ? "#define USE_MORPHNORMALS" : "", i.doubleSided ? "#define DOUBLE_SIDED" : "", i.flipSided ? "#define FLIP_SIDED" : "", i.shadowMapEnabled ? "#define USE_SHADOWMAP" : "", i.shadowMapEnabled ? "#define " + l : "", i.sizeAttenuation ? "#define USE_SIZEATTENUATION" : "", i.logarithmicDepthBuffer ? "#define USE_LOGDEPTHBUF" : "", i.logarithmicDepthBuffer && i.rendererExtensionFragDepth ? "#define USE_LOGDEPTHBUF_EXT" : "", "uniform mat4 modelMatrix;", "uniform mat4 modelViewMatrix;", "uniform mat4 projectionMatrix;", "uniform mat4 viewMatrix;", "uniform mat3 normalMatrix;", "uniform vec3 cameraPosition;", "uniform bool isOrthographic;", "#ifdef USE_INSTANCING", "\tattribute mat4 instanceMatrix;", "#endif", "#ifdef USE_INSTANCING_COLOR", "\tattribute vec3 instanceColor;", "#endif", "attribute vec3 position;", "attribute vec3 normal;", "attribute vec2 uv;", "#ifdef USE_TANGENT", "\tattribute vec4 tangent;", "#endif", "#ifdef USE_COLOR", "\tattribute vec3 color;", "#endif", "#ifdef USE_MORPHTARGETS", "\tattribute vec3 morphTarget0;", "\tattribute vec3 morphTarget1;", "\tattribute vec3 morphTarget2;", "\tattribute vec3 morphTarget3;", "\t#ifdef USE_MORPHNORMALS", "\t\tattribute vec3 morphNormal0;", "\t\tattribute vec3 morphNormal1;", "\t\tattribute vec3 morphNormal2;", "\t\tattribute vec3 morphNormal3;", "\t#else", "\t\tattribute vec3 morphTarget4;", "\t\tattribute vec3 morphTarget5;", "\t\tattribute vec3 morphTarget6;", "\t\tattribute vec3 morphTarget7;", "\t#endif", "#endif", "#ifdef USE_SKINNING", "\tattribute vec4 skinIndex;", "\tattribute vec4 skinWeight;", "#endif", "\n"].filter(filterEmptyLine).join("\n"), f = [p, generatePrecision(i), "#define SHADER_NAME " + i.shaderName, m, i.alphaTest ? "#define ALPHATEST " + i.alphaTest + (i.alphaTest % 1 ? "" : ".0") : "", "#define GAMMA_FACTOR " + d, i.useFog && i.fog ? "#define USE_FOG" : "", i.useFog && i.fogExp2 ? "#define FOG_EXP2" : "", i.map ? "#define USE_MAP" : "", i.matcap ? "#define USE_MATCAP" : "", i.envMap ? "#define USE_ENVMAP" : "", i.envMap ? "#define " + c : "", i.envMap ? "#define " + h : "", i.envMap ? "#define " + u : "", i.lightMap ? "#define USE_LIGHTMAP" : "", i.aoMap ? "#define USE_AOMAP" : "", i.emissiveMap ? "#define USE_EMISSIVEMAP" : "", i.bumpMap ? "#define USE_BUMPMAP" : "", i.normalMap ? "#define USE_NORMALMAP" : "", i.normalMap && i.objectSpaceNormalMap ? "#define OBJECTSPACE_NORMALMAP" : "", i.normalMap && i.tangentSpaceNormalMap ? "#define TANGENTSPACE_NORMALMAP" : "", i.clearcoatMap ? "#define USE_CLEARCOATMAP" : "", i.clearcoatRoughnessMap ? "#define USE_CLEARCOAT_ROUGHNESSMAP" : "", i.clearcoatNormalMap ? "#define USE_CLEARCOAT_NORMALMAP" : "", i.specularMap ? "#define USE_SPECULARMAP" : "", i.roughnessMap ? "#define USE_ROUGHNESSMAP" : "", i.metalnessMap ? "#define USE_METALNESSMAP" : "", i.alphaMap ? "#define USE_ALPHAMAP" : "", i.sheen ? "#define USE_SHEEN" : "", i.transmissionMap ? "#define USE_TRANSMISSIONMAP" : "", i.vertexTangents ? "#define USE_TANGENT" : "", i.vertexColors || i.instancingColor ? "#define USE_COLOR" : "", i.vertexUvs ? "#define USE_UV" : "", i.uvsVertexOnly ? "#define UVS_VERTEX_ONLY" : "", i.gradientMap ? "#define USE_GRADIENTMAP" : "", i.flatShading ? "#define FLAT_SHADED" : "", i.doubleSided ? "#define DOUBLE_SIDED" : "", i.flipSided ? "#define FLIP_SIDED" : "", i.shadowMapEnabled ? "#define USE_SHADOWMAP" : "", i.shadowMapEnabled ? "#define " + l : "", i.premultipliedAlpha ? "#define PREMULTIPLIED_ALPHA" : "", i.physicallyCorrectLights ? "#define PHYSICALLY_CORRECT_LIGHTS" : "", i.logarithmicDepthBuffer ? "#define USE_LOGDEPTHBUF" : "", i.logarithmicDepthBuffer && i.rendererExtensionFragDepth ? "#define USE_LOGDEPTHBUF_EXT" : "", (i.extensionShaderTextureLOD || i.envMap) && i.rendererExtensionShaderTextureLod ? "#define TEXTURE_LOD_EXT" : "", "uniform mat4 viewMatrix;", "uniform vec3 cameraPosition;", "uniform bool isOrthographic;", 0 !== i.toneMapping ? "#define TONE_MAPPING" : "", 0 !== i.toneMapping ? ShaderChunk.tonemapping_pars_fragment : "", 0 !== i.toneMapping ? getToneMappingFunction("toneMapping", i.toneMapping) : "", i.dithering ? "#define DITHERING" : "", ShaderChunk.encodings_pars_fragment, i.map ? getTexelDecodingFunction("mapTexelToLinear", i.mapEncoding) : "", i.matcap ? getTexelDecodingFunction("matcapTexelToLinear", i.matcapEncoding) : "", i.envMap ? getTexelDecodingFunction("envMapTexelToLinear", i.envMapEncoding) : "", i.emissiveMap ? getTexelDecodingFunction("emissiveMapTexelToLinear", i.emissiveMapEncoding) : "", i.lightMap ? getTexelDecodingFunction("lightMapTexelToLinear", i.lightMapEncoding) : "", getTexelEncodingFunction("linearToOutputTexel", i.outputEncoding), i.depthPacking ? "#define DEPTH_PACKING " + i.depthPacking : "", "\n"].filter(filterEmptyLine).join("\n")), s = resolveIncludes(s), s = replaceLightNums(s, i), s = replaceClippingPlaneNums(s, i), o = resolveIncludes(o), o = replaceLightNums(o, i), o = replaceClippingPlaneNums(o, i), s = unrollLoops(s), o = unrollLoops(o), i.isWebGL2 && !0 !== i.isRawShaderMaterial && (v = "#version 300 es\n", g = ["#define attribute in", "#define varying out", "#define texture2D texture"].join("\n") + "\n" + g, f = ["#define varying in", i.glslVersion === GLSL3 ? "" : "out highp vec4 pc_fragColor;", i.glslVersion === GLSL3 ? "" : "#define gl_FragColor pc_fragColor", "#define gl_FragDepthEXT gl_FragDepth", "#define texture2D texture", "#define textureCube texture", "#define texture2DProj textureProj", "#define texture2DLodEXT textureLod", "#define texture2DProjLodEXT textureProjLod", "#define textureCubeLodEXT textureLod", "#define texture2DGradEXT textureGrad", "#define texture2DProjGradEXT textureProjGrad", "#define textureCubeGradEXT textureGrad"].join("\n") + "\n" + f);
	const y = v + f + o,
		E = WebGLShader(r, 35633, v + g + s),
		_ = WebGLShader(r, 35632, y);
	if (r.attachShader(A, E), r.attachShader(A, _), void 0 !== i.index0AttributeName ? r.bindAttribLocation(A, 0, i.index0AttributeName) : !0 === i.morphTargets && r.bindAttribLocation(A, 0, "position"), r.linkProgram(A), e.debug.checkShaderErrors) {
		const e = r.getProgramInfoLog(A).trim(),
			t = r.getShaderInfoLog(E).trim(),
			i = r.getShaderInfoLog(_).trim();
		let n = !0,
			a = !0;
		if (!1 === r.getProgramParameter(A, 35714)) {
			n = !1;
			const t = getShaderErrors(r, E, "vertex"),
				i = getShaderErrors(r, _, "fragment");
			console.error("THREE.WebGLProgram: shader error: ", r.getError(), "35715", r.getProgramParameter(A, 35715), "gl.getProgramInfoLog", e, t, i)
		} else "" !== e ? console.warn("THREE.WebGLProgram: gl.getProgramInfoLog()", e) : "" !== t && "" !== i || (a = !1);
		a && (this.diagnostics = {
			runnable: n,
			programLog: e,
			vertexShader: {
				log: t,
				prefix: g
			},
			fragmentShader: {
				log: i,
				prefix: f
			}
		})
	}
	let b, x;
	return r.deleteShader(E), r.deleteShader(_), this.getUniforms = function() {
		return void 0 === b && (b = new WebGLUniforms(r, A)), b
	}, this.getAttributes = function() {
		return void 0 === x && (x = fetchAttributeLocations(r, A)), x
	}, this.destroy = function() {
		n.releaseStatesOfProgram(this), r.deleteProgram(A), this.program = void 0
	}, this.name = i.shaderName, this.id = programIdCount++, this.cacheKey = t, this.usedTimes = 1, this.program = A, this.vertexShader = E, this.fragmentShader = _, this
}

function WebGLPrograms(e, t, i, n, r, a) {
	const s = [],
		o = n.isWebGL2,
		l = n.logarithmicDepthBuffer,
		c = n.floatVertexTextures,
		h = n.maxVertexUniforms,
		u = n.vertexTextures;
	let d = n.precision;
	const p = {
			MeshDepthMaterial: "depth",
			MeshDistanceMaterial: "distanceRGBA",
			MeshNormalMaterial: "normal",
			MeshBasicMaterial: "basic",
			MeshLambertMaterial: "lambert",
			MeshPhongMaterial: "phong",
			MeshToonMaterial: "toon",
			MeshStandardMaterial: "physical",
			MeshPhysicalMaterial: "physical",
			MeshMatcapMaterial: "matcap",
			LineBasicMaterial: "basic",
			LineDashedMaterial: "dashed",
			PointsMaterial: "points",
			ShadowMaterial: "shadow",
			SpriteMaterial: "sprite"
		},
		m = ["precision", "isWebGL2", "supportsVertexTextures", "outputEncoding", "instancing", "instancingColor", "map", "mapEncoding", "matcap", "matcapEncoding", "envMap", "envMapMode", "envMapEncoding", "envMapCubeUV", "lightMap", "lightMapEncoding", "aoMap", "emissiveMap", "emissiveMapEncoding", "bumpMap", "normalMap", "objectSpaceNormalMap", "tangentSpaceNormalMap", "clearcoatMap", "clearcoatRoughnessMap", "clearcoatNormalMap", "displacementMap", "specularMap", "roughnessMap", "metalnessMap", "gradientMap", "alphaMap", "combine", "vertexColors", "vertexTangents", "vertexUvs", "uvsVertexOnly", "fog", "useFog", "fogExp2", "flatShading", "sizeAttenuation", "logarithmicDepthBuffer", "skinning", "maxBones", "useVertexTexture", "morphTargets", "morphNormals", "maxMorphTargets", "maxMorphNormals", "premultipliedAlpha", "numDirLights", "numPointLights", "numSpotLights", "numHemiLights", "numRectAreaLights", "numDirLightShadows", "numPointLightShadows", "numSpotLightShadows", "shadowMapEnabled", "shadowMapType", "toneMapping", "physicallyCorrectLights", "alphaTest", "doubleSided", "flipSided", "numClippingPlanes", "numClipIntersection", "depthPacking", "dithering", "sheen", "transmissionMap"];

	function A(e) {
		let t;
		return e && e.isTexture ? t = e.encoding : e && e.isWebGLRenderTarget ? (console.warn("THREE.WebGLPrograms.getTextureEncodingFromMap: don't use render targets as textures. Use their .texture property instead."), t = e.texture.encoding) : t = 3e3, t
	}
	return {
		getParameters: function(r, s, m, g, f) {
			const v = g.fog,
				y = r.isMeshStandardMaterial ? g.environment : null,
				E = t.get(r.envMap || y),
				_ = p[r.type],
				b = f.isSkinnedMesh ? function(e) {
					const t = e.skeleton.bones;
					if (c) return 1024; {
						const e = h,
							i = Math.floor((e - 20) / 4),
							n = Math.min(i, t.length);
						return n < t.length ? (console.warn("THREE.WebGLRenderer: Skeleton has " + t.length + " bones. This GPU supports " + n + "."), 0) : n
					}
				}(f) : 0;
			let x, w;
			if (null !== r.precision && (d = n.getMaxPrecision(r.precision), d !== r.precision && console.warn("THREE.WebGLProgram.getParameters:", r.precision, "not supported, using", d, "instead.")), _) {
				const e = ShaderLib[_];
				x = e.vertexShader, w = e.fragmentShader
			} else x = r.vertexShader, w = r.fragmentShader;
			const C = e.getRenderTarget();
			return {
				isWebGL2: o,
				shaderID: _,
				shaderName: r.type,
				vertexShader: x,
				fragmentShader: w,
				defines: r.defines,
				isRawShaderMaterial: !0 === r.isRawShaderMaterial,
				glslVersion: r.glslVersion,
				precision: d,
				instancing: !0 === f.isInstancedMesh,
				instancingColor: !0 === f.isInstancedMesh && null !== f.instanceColor,
				supportsVertexTextures: u,
				outputEncoding: null !== C ? A(C.texture) : e.outputEncoding,
				map: !!r.map,
				mapEncoding: A(r.map),
				matcap: !!r.matcap,
				matcapEncoding: A(r.matcap),
				envMap: !!E,
				envMapMode: E && E.mapping,
				envMapEncoding: A(E),
				envMapCubeUV: !!E && (306 === E.mapping || 307 === E.mapping),
				lightMap: !!r.lightMap,
				lightMapEncoding: A(r.lightMap),
				aoMap: !!r.aoMap,
				emissiveMap: !!r.emissiveMap,
				emissiveMapEncoding: A(r.emissiveMap),
				bumpMap: !!r.bumpMap,
				normalMap: !!r.normalMap,
				objectSpaceNormalMap: 1 === r.normalMapType,
				tangentSpaceNormalMap: 0 === r.normalMapType,
				clearcoatMap: !!r.clearcoatMap,
				clearcoatRoughnessMap: !!r.clearcoatRoughnessMap,
				clearcoatNormalMap: !!r.clearcoatNormalMap,
				displacementMap: !!r.displacementMap,
				roughnessMap: !!r.roughnessMap,
				metalnessMap: !!r.metalnessMap,
				specularMap: !!r.specularMap,
				alphaMap: !!r.alphaMap,
				gradientMap: !!r.gradientMap,
				sheen: !!r.sheen,
				transmissionMap: !!r.transmissionMap,
				combine: r.combine,
				vertexTangents: r.normalMap && r.vertexTangents,
				vertexColors: r.vertexColors,
				vertexUvs: !!(r.map || r.bumpMap || r.normalMap || r.specularMap || r.alphaMap || r.emissiveMap || r.roughnessMap || r.metalnessMap || r.clearcoatMap || r.clearcoatRoughnessMap || r.clearcoatNormalMap || r.displacementMap || r.transmissionMap),
				uvsVertexOnly: !(r.map || r.bumpMap || r.normalMap || r.specularMap || r.alphaMap || r.emissiveMap || r.roughnessMap || r.metalnessMap || r.clearcoatNormalMap || r.transmissionMap || !r.displacementMap),
				fog: !!v,
				useFog: r.fog,
				fogExp2: v && v.isFogExp2,
				flatShading: !!r.flatShading,
				sizeAttenuation: r.sizeAttenuation,
				logarithmicDepthBuffer: l,
				skinning: r.skinning && b > 0,
				maxBones: b,
				useVertexTexture: c,
				morphTargets: r.morphTargets,
				morphNormals: r.morphNormals,
				maxMorphTargets: e.maxMorphTargets,
				maxMorphNormals: e.maxMorphNormals,
				numDirLights: s.directional.length,
				numPointLights: s.point.length,
				numSpotLights: s.spot.length,
				numRectAreaLights: s.rectArea.length,
				numHemiLights: s.hemi.length,
				numDirLightShadows: s.directionalShadowMap.length,
				numPointLightShadows: s.pointShadowMap.length,
				numSpotLightShadows: s.spotShadowMap.length,
				numClippingPlanes: a.numPlanes,
				numClipIntersection: a.numIntersection,
				dithering: r.dithering,
				shadowMapEnabled: e.shadowMap.enabled && m.length > 0,
				shadowMapType: e.shadowMap.type,
				toneMapping: r.toneMapped ? e.toneMapping : 0,
				physicallyCorrectLights: e.physicallyCorrectLights,
				premultipliedAlpha: r.premultipliedAlpha,
				alphaTest: r.alphaTest,
				doubleSided: 2 === r.side,
				flipSided: 1 === r.side,
				depthPacking: void 0 !== r.depthPacking && r.depthPacking,
				index0AttributeName: r.index0AttributeName,
				extensionDerivatives: r.extensions && r.extensions.derivatives,
				extensionFragDepth: r.extensions && r.extensions.fragDepth,
				extensionDrawBuffers: r.extensions && r.extensions.drawBuffers,
				extensionShaderTextureLOD: r.extensions && r.extensions.shaderTextureLOD,
				rendererExtensionFragDepth: o || i.has("EXT_frag_depth"),
				rendererExtensionDrawBuffers: o || i.has("WEBGL_draw_buffers"),
				rendererExtensionShaderTextureLod: o || i.has("EXT_shader_texture_lod"),
				customProgramCacheKey: r.customProgramCacheKey()
			}
		},
		getProgramCacheKey: function(t) {
			const i = [];
			if (t.shaderID ? i.push(t.shaderID) : (i.push(t.fragmentShader), i.push(t.vertexShader)), void 0 !== t.defines)
				for (const e in t.defines) i.push(e), i.push(t.defines[e]);
			if (!1 === t.isRawShaderMaterial) {
				for (let e = 0; e < m.length; e++) i.push(t[m[e]]);
				i.push(e.outputEncoding), i.push(e.gammaFactor)
			}
			return i.push(t.customProgramCacheKey), i.join()
		},
		getUniforms: function(e) {
			const t = p[e.type];
			let i;
			if (t) {
				const e = ShaderLib[t];
				i = UniformsUtils.clone(e.uniforms)
			} else i = e.uniforms;
			return i
		},
		acquireProgram: function(t, i) {
			let n;
			for (let e = 0, t = s.length; e < t; e++) {
				const t = s[e];
				if (t.cacheKey === i) {
					n = t, ++n.usedTimes;
					break
				}
			}
			return void 0 === n && (n = new WebGLProgram(e, i, t, r), s.push(n)), n
		},
		releaseProgram: function(e) {
			if (0 == --e.usedTimes) {
				const t = s.indexOf(e);
				s[t] = s[s.length - 1], s.pop(), e.destroy()
			}
		},
		programs: s
	}
}

function WebGLProperties() {
	let e = new WeakMap;
	return {
		get: function(t) {
			let i = e.get(t);
			return void 0 === i && (i = {}, e.set(t, i)), i
		},
		remove: function(t) {
			e.delete(t)
		},
		update: function(t, i, n) {
			e.get(t)[i] = n
		},
		dispose: function() {
			e = new WeakMap
		}
	}
}

function painterSortStable(e, t) {
	return e.groupOrder !== t.groupOrder ? e.groupOrder - t.groupOrder : e.renderOrder !== t.renderOrder ? e.renderOrder - t.renderOrder : e.program !== t.program ? e.program.id - t.program.id : e.material.id !== t.material.id ? e.material.id - t.material.id : e.z !== t.z ? e.z - t.z : e.id - t.id
}

function reversePainterSortStable(e, t) {
	return e.groupOrder !== t.groupOrder ? e.groupOrder - t.groupOrder : e.renderOrder !== t.renderOrder ? e.renderOrder - t.renderOrder : e.z !== t.z ? t.z - e.z : e.id - t.id
}

function WebGLRenderList(e) {
	const t = [];
	let i = 0;
	const n = [],
		r = [],
		a = {
			id: -1
		};

	function s(n, r, s, o, l, c) {
		let h = t[i];
		const u = e.get(s);
		return void 0 === h ? (h = {
			id: n.id,
			object: n,
			geometry: r,
			material: s,
			program: u.program || a,
			groupOrder: o,
			renderOrder: n.renderOrder,
			z: l,
			group: c
		}, t[i] = h) : (h.id = n.id, h.object = n, h.geometry = r, h.material = s, h.program = u.program || a, h.groupOrder = o, h.renderOrder = n.renderOrder, h.z = l, h.group = c), i++, h
	}
	return {
		opaque: n,
		transparent: r,
		init: function() {
			i = 0, n.length = 0, r.length = 0
		},
		push: function(e, t, i, a, o, l) {
			const c = s(e, t, i, a, o, l);
			(!0 === i.transparent ? r : n).push(c)
		},
		unshift: function(e, t, i, a, o, l) {
			const c = s(e, t, i, a, o, l);
			(!0 === i.transparent ? r : n).unshift(c)
		},
		finish: function() {
			for (let e = i, n = t.length; e < n; e++) {
				const i = t[e];
				if (null === i.id) break;
				i.id = null, i.object = null, i.geometry = null, i.material = null, i.program = null, i.group = null
			}
		},
		sort: function(e, t) {
			n.length > 1 && n.sort(e || painterSortStable), r.length > 1 && r.sort(t || reversePainterSortStable)
		}
	}
}

function WebGLRenderLists(e) {
	let t = new WeakMap;
	return {
		get: function(i, n) {
			let r;
			return !1 === t.has(i) ? (r = new WebGLRenderList(e), t.set(i, [r])) : n >= t.get(i).length ? (r = new WebGLRenderList(e), t.get(i).push(r)) : r = t.get(i)[n], r
		},
		dispose: function() {
			t = new WeakMap
		}
	}
}

function UniformsCache() {
	const e = {};
	return {
		get: function(t) {
			if (void 0 !== e[t.id]) return e[t.id];
			let i;
			switch (t.type) {
				case "DirectionalLight":
					i = {
						direction: new Vector3,
						color: new Color
					};
					break;
				case "SpotLight":
					i = {
						position: new Vector3,
						direction: new Vector3,
						color: new Color,
						distance: 0,
						coneCos: 0,
						penumbraCos: 0,
						decay: 0
					};
					break;
				case "PointLight":
					i = {
						position: new Vector3,
						color: new Color,
						distance: 0,
						decay: 0
					};
					break;
				case "HemisphereLight":
					i = {
						direction: new Vector3,
						skyColor: new Color,
						groundColor: new Color
					};
					break;
				case "RectAreaLight":
					i = {
						color: new Color,
						position: new Vector3,
						halfWidth: new Vector3,
						halfHeight: new Vector3
					}
			}
			return e[t.id] = i, i
		}
	}
}

function ShadowUniformsCache() {
	const e = {};
	return {
		get: function(t) {
			if (void 0 !== e[t.id]) return e[t.id];
			let i;
			switch (t.type) {
				case "DirectionalLight":
				case "SpotLight":
					i = {
						shadowBias: 0,
						shadowNormalBias: 0,
						shadowRadius: 1,
						shadowMapSize: new Vector2
					};
					break;
				case "PointLight":
					i = {
						shadowBias: 0,
						shadowNormalBias: 0,
						shadowRadius: 1,
						shadowMapSize: new Vector2,
						shadowCameraNear: 1,
						shadowCameraFar: 1e3
					}
			}
			return e[t.id] = i, i
		}
	}
}
let nextVersion = 0;

function shadowCastingLightsFirst(e, t) {
	return (t.castShadow ? 1 : 0) - (e.castShadow ? 1 : 0)
}

function WebGLLights(e, t) {
	const i = new UniformsCache,
		n = ShadowUniformsCache(),
		r = {
			version: 0,
			hash: {
				directionalLength: -1,
				pointLength: -1,
				spotLength: -1,
				rectAreaLength: -1,
				hemiLength: -1,
				numDirectionalShadows: -1,
				numPointShadows: -1,
				numSpotShadows: -1
			},
			ambient: [0, 0, 0],
			probe: [],
			directional: [],
			directionalShadow: [],
			directionalShadowMap: [],
			directionalShadowMatrix: [],
			spot: [],
			spotShadow: [],
			spotShadowMap: [],
			spotShadowMatrix: [],
			rectArea: [],
			rectAreaLTC1: null,
			rectAreaLTC2: null,
			point: [],
			pointShadow: [],
			pointShadowMap: [],
			pointShadowMatrix: [],
			hemi: []
		};
	for (let e = 0; e < 9; e++) r.probe.push(new Vector3);
	const a = new Vector3,
		s = new Matrix4,
		o = new Matrix4;
	return {
		setup: function(a) {
			let s = 0,
				o = 0,
				l = 0;
			for (let e = 0; e < 9; e++) r.probe[e].set(0, 0, 0);
			let c = 0,
				h = 0,
				u = 0,
				d = 0,
				p = 0,
				m = 0,
				A = 0,
				g = 0;
			a.sort(shadowCastingLightsFirst);
			for (let e = 0, t = a.length; e < t; e++) {
				const t = a[e],
					f = t.color,
					v = t.intensity,
					y = t.distance,
					E = t.shadow && t.shadow.map ? t.shadow.map.texture : null;
				if (t.isAmbientLight) s += f.r * v, o += f.g * v, l += f.b * v;
				else if (t.isLightProbe)
					for (let e = 0; e < 9; e++) r.probe[e].addScaledVector(t.sh.coefficients[e], v);
				else if (t.isDirectionalLight) {
					const e = i.get(t);
					if (e.color.copy(t.color).multiplyScalar(t.intensity), t.castShadow) {
						const e = t.shadow,
							i = n.get(t);
						i.shadowBias = e.bias, i.shadowNormalBias = e.normalBias, i.shadowRadius = e.radius, i.shadowMapSize = e.mapSize, r.directionalShadow[c] = i, r.directionalShadowMap[c] = E, r.directionalShadowMatrix[c] = t.shadow.matrix, m++
					}
					r.directional[c] = e, c++
				} else if (t.isSpotLight) {
					const e = i.get(t);
					if (e.position.setFromMatrixPosition(t.matrixWorld), e.color.copy(f).multiplyScalar(v), e.distance = y, e.coneCos = Math.cos(t.angle), e.penumbraCos = Math.cos(t.angle * (1 - t.penumbra)), e.decay = t.decay, t.castShadow) {
						const e = t.shadow,
							i = n.get(t);
						i.shadowBias = e.bias, i.shadowNormalBias = e.normalBias, i.shadowRadius = e.radius, i.shadowMapSize = e.mapSize, r.spotShadow[u] = i, r.spotShadowMap[u] = E, r.spotShadowMatrix[u] = t.shadow.matrix, g++
					}
					r.spot[u] = e, u++
				} else if (t.isRectAreaLight) {
					const e = i.get(t);
					e.color.copy(f).multiplyScalar(v), e.halfWidth.set(.5 * t.width, 0, 0), e.halfHeight.set(0, .5 * t.height, 0), r.rectArea[d] = e, d++
				} else if (t.isPointLight) {
					const e = i.get(t);
					if (e.color.copy(t.color).multiplyScalar(t.intensity), e.distance = t.distance, e.decay = t.decay, t.castShadow) {
						const e = t.shadow,
							i = n.get(t);
						i.shadowBias = e.bias, i.shadowNormalBias = e.normalBias, i.shadowRadius = e.radius, i.shadowMapSize = e.mapSize, i.shadowCameraNear = e.camera.near, i.shadowCameraFar = e.camera.far, r.pointShadow[h] = i, r.pointShadowMap[h] = E, r.pointShadowMatrix[h] = t.shadow.matrix, A++
					}
					r.point[h] = e, h++
				} else if (t.isHemisphereLight) {
					const e = i.get(t);
					e.skyColor.copy(t.color).multiplyScalar(v), e.groundColor.copy(t.groundColor).multiplyScalar(v), r.hemi[p] = e, p++
				}
			}
			d > 0 && (t.isWebGL2 || !0 === e.has("OES_texture_float_linear") ? (r.rectAreaLTC1 = UniformsLib.LTC_FLOAT_1, r.rectAreaLTC2 = UniformsLib.LTC_FLOAT_2) : !0 === e.has("OES_texture_half_float_linear") ? (r.rectAreaLTC1 = UniformsLib.LTC_HALF_1, r.rectAreaLTC2 = UniformsLib.LTC_HALF_2) : console.error("THREE.WebGLRenderer: Unable to use RectAreaLight. Missing WebGL extensions.")), r.ambient[0] = s, r.ambient[1] = o, r.ambient[2] = l;
			const f = r.hash;
			f.directionalLength === c && f.pointLength === h && f.spotLength === u && f.rectAreaLength === d && f.hemiLength === p && f.numDirectionalShadows === m && f.numPointShadows === A && f.numSpotShadows === g || (r.directional.length = c, r.spot.length = u, r.rectArea.length = d, r.point.length = h, r.hemi.length = p, r.directionalShadow.length = m, r.directionalShadowMap.length = m, r.pointShadow.length = A, r.pointShadowMap.length = A, r.spotShadow.length = g, r.spotShadowMap.length = g, r.directionalShadowMatrix.length = m, r.pointShadowMatrix.length = A, r.spotShadowMatrix.length = g, f.directionalLength = c, f.pointLength = h, f.spotLength = u, f.rectAreaLength = d, f.hemiLength = p, f.numDirectionalShadows = m, f.numPointShadows = A, f.numSpotShadows = g, r.version = nextVersion++)
		},
		setupView: function(e, t) {
			let i = 0,
				n = 0,
				l = 0,
				c = 0,
				h = 0;
			const u = t.matrixWorldInverse;
			for (let t = 0, d = e.length; t < d; t++) {
				const d = e[t];
				if (d.isDirectionalLight) {
					const e = r.directional[i];
					e.direction.setFromMatrixPosition(d.matrixWorld), a.setFromMatrixPosition(d.target.matrixWorld), e.direction.sub(a), e.direction.transformDirection(u), i++
				} else if (d.isSpotLight) {
					const e = r.spot[l];
					e.position.setFromMatrixPosition(d.matrixWorld), e.position.applyMatrix4(u), e.direction.setFromMatrixPosition(d.matrixWorld), a.setFromMatrixPosition(d.target.matrixWorld), e.direction.sub(a), e.direction.transformDirection(u), l++
				} else if (d.isRectAreaLight) {
					const e = r.rectArea[c];
					e.position.setFromMatrixPosition(d.matrixWorld), e.position.applyMatrix4(u), o.identity(), s.copy(d.matrixWorld), s.premultiply(u), o.extractRotation(s), e.halfWidth.set(.5 * d.width, 0, 0), e.halfHeight.set(0, .5 * d.height, 0), e.halfWidth.applyMatrix4(o), e.halfHeight.applyMatrix4(o), c++
				} else if (d.isPointLight) {
					const e = r.point[n];
					e.position.setFromMatrixPosition(d.matrixWorld), e.position.applyMatrix4(u), n++
				} else if (d.isHemisphereLight) {
					const e = r.hemi[h];
					e.direction.setFromMatrixPosition(d.matrixWorld), e.direction.transformDirection(u), e.direction.normalize(), h++
				}
			}
		},
		state: r
	}
}

function WebGLRenderState(e, t) {
	const i = new WebGLLights(e, t),
		n = [],
		r = [];
	return {
		init: function() {
			n.length = 0, r.length = 0
		},
		state: {
			lightsArray: n,
			shadowsArray: r,
			lights: i
		},
		setupLights: function() {
			i.setup(n)
		},
		setupLightsView: function(e) {
			i.setupView(n, e)
		},
		pushLight: function(e) {
			n.push(e)
		},
		pushShadow: function(e) {
			r.push(e)
		}
	}
}

function WebGLRenderStates(e, t) {
	let i = new WeakMap;
	return {
		get: function(n, r = 0) {
			let a;
			return !1 === i.has(n) ? (a = new WebGLRenderState(e, t), i.set(n, [a])) : r >= i.get(n).length ? (a = new WebGLRenderState(e, t), i.get(n).push(a)) : a = i.get(n)[r], a
		},
		dispose: function() {
			i = new WeakMap
		}
	}
}
class MeshDepthMaterial extends Material$1 {
	constructor(e) {
		super(), this.type = "MeshDepthMaterial", this.depthPacking = 3200, this.skinning = !1, this.morphTargets = !1, this.map = null, this.alphaMap = null, this.displacementMap = null, this.displacementScale = 1, this.displacementBias = 0, this.wireframe = !1, this.wireframeLinewidth = 1, this.fog = !1, this.setValues(e)
	}
	copy(e) {
		return super.copy(e), this.depthPacking = e.depthPacking, this.skinning = e.skinning, this.morphTargets = e.morphTargets, this.map = e.map, this.alphaMap = e.alphaMap, this.displacementMap = e.displacementMap, this.displacementScale = e.displacementScale, this.displacementBias = e.displacementBias, this.wireframe = e.wireframe, this.wireframeLinewidth = e.wireframeLinewidth, this
	}
}
MeshDepthMaterial.prototype.isMeshDepthMaterial = !0;
class MeshDistanceMaterial extends Material$1 {
	constructor(e) {
		super(), this.type = "MeshDistanceMaterial", this.referencePosition = new Vector3, this.nearDistance = 1, this.farDistance = 1e3, this.skinning = !1, this.morphTargets = !1, this.map = null, this.alphaMap = null, this.displacementMap = null, this.displacementScale = 1, this.displacementBias = 0, this.fog = !1, this.setValues(e)
	}
	copy(e) {
		return super.copy(e), this.referencePosition.copy(e.referencePosition), this.nearDistance = e.nearDistance, this.farDistance = e.farDistance, this.skinning = e.skinning, this.morphTargets = e.morphTargets, this.map = e.map, this.alphaMap = e.alphaMap, this.displacementMap = e.displacementMap, this.displacementScale = e.displacementScale, this.displacementBias = e.displacementBias, this
	}
}
MeshDistanceMaterial.prototype.isMeshDistanceMaterial = !0;
var vsm_frag = "uniform sampler2D shadow_pass;\nuniform vec2 resolution;\nuniform float radius;\n#include <packing>\nvoid main() {\n\tfloat mean = 0.0;\n\tfloat squared_mean = 0.0;\n\tfloat depth = unpackRGBAToDepth( texture2D( shadow_pass, ( gl_FragCoord.xy ) / resolution ) );\n\tfor ( float i = -1.0; i < 1.0 ; i += SAMPLE_RATE) {\n\t\t#ifdef HORIZONTAL_PASS\n\t\t\tvec2 distribution = unpackRGBATo2Half( texture2D( shadow_pass, ( gl_FragCoord.xy + vec2( i, 0.0 ) * radius ) / resolution ) );\n\t\t\tmean += distribution.x;\n\t\t\tsquared_mean += distribution.y * distribution.y + distribution.x * distribution.x;\n\t\t#else\n\t\t\tfloat depth = unpackRGBAToDepth( texture2D( shadow_pass, ( gl_FragCoord.xy + vec2( 0.0, i ) * radius ) / resolution ) );\n\t\t\tmean += depth;\n\t\t\tsquared_mean += depth * depth;\n\t\t#endif\n\t}\n\tmean = mean * HALF_SAMPLE_RATE;\n\tsquared_mean = squared_mean * HALF_SAMPLE_RATE;\n\tfloat std_dev = sqrt( squared_mean - mean * mean );\n\tgl_FragColor = pack2HalfToRGBA( vec2( mean, std_dev ) );\n}",
	vsm_vert = "void main() {\n\tgl_Position = vec4( position, 1.0 );\n}";

function WebGLShadowMap(e, t, i) {
	let n = new Frustum;
	const r = new Vector2,
		a = new Vector2,
		s = new Vector4,
		o = [],
		l = [],
		c = {},
		h = {
			0: 1,
			1: 0,
			2: 2
		},
		u = new ShaderMaterial({
			defines: {
				SAMPLE_RATE: 2 / 8,
				HALF_SAMPLE_RATE: 1 / 8
			},
			uniforms: {
				shadow_pass: {
					value: null
				},
				resolution: {
					value: new Vector2
				},
				radius: {
					value: 4
				}
			},
			vertexShader: vsm_vert,
			fragmentShader: vsm_frag
		}),
		d = u.clone();
	d.defines.HORIZONTAL_PASS = 1;
	const p = new BufferGeometry;
	p.setAttribute("position", new BufferAttribute(new Float32Array([-1, -1, .5, 3, -1, .5, -1, 3, .5]), 3));
	const m = new Mesh(p, u),
		A = this;

	function g(i, n) {
		const r = t.update(m);
		u.uniforms.shadow_pass.value = i.map.texture, u.uniforms.resolution.value = i.mapSize, u.uniforms.radius.value = i.radius, e.setRenderTarget(i.mapPass), e.clear(), e.renderBufferDirect(n, null, r, u, m, null), d.uniforms.shadow_pass.value = i.mapPass.texture, d.uniforms.resolution.value = i.mapSize, d.uniforms.radius.value = i.radius, e.setRenderTarget(i.map), e.clear(), e.renderBufferDirect(n, null, r, d, m, null)
	}

	function f(e, t, i) {
		const n = e << 0 | t << 1 | i << 2;
		let r = o[n];
		return void 0 === r && (r = new MeshDepthMaterial({
			depthPacking: 3201,
			morphTargets: e,
			skinning: t
		}), o[n] = r), r
	}

	function v(e, t, i) {
		const n = e << 0 | t << 1 | i << 2;
		let r = l[n];
		return void 0 === r && (r = new MeshDistanceMaterial({
			morphTargets: e,
			skinning: t
		}), l[n] = r), r
	}

	function y(t, i, n, r, a, s, o) {
		let l = null,
			u = f,
			d = t.customDepthMaterial;
		if (!0 === r.isPointLight && (u = v, d = t.customDistanceMaterial), void 0 === d) {
			let e = !1;
			!0 === n.morphTargets && (e = i.morphAttributes && i.morphAttributes.position && i.morphAttributes.position.length > 0);
			let r = !1;
			!0 === t.isSkinnedMesh && (!0 === n.skinning ? r = !0 : console.warn("THREE.WebGLShadowMap: THREE.SkinnedMesh with material.skinning set to false:", t));
			l = u(e, r, !0 === t.isInstancedMesh)
		} else l = d;
		if (e.localClippingEnabled && !0 === n.clipShadows && 0 !== n.clippingPlanes.length) {
			const e = l.uuid,
				t = n.uuid;
			let i = c[e];
			void 0 === i && (i = {}, c[e] = i);
			let r = i[t];
			void 0 === r && (r = l.clone(), i[t] = r), l = r
		}
		return l.visible = n.visible, l.wireframe = n.wireframe, l.side = 3 === o ? null !== n.shadowSide ? n.shadowSide : n.side : null !== n.shadowSide ? n.shadowSide : h[n.side], l.clipShadows = n.clipShadows, l.clippingPlanes = n.clippingPlanes, l.clipIntersection = n.clipIntersection, l.wireframeLinewidth = n.wireframeLinewidth, l.linewidth = n.linewidth, !0 === r.isPointLight && !0 === l.isMeshDistanceMaterial && (l.referencePosition.setFromMatrixPosition(r.matrixWorld), l.nearDistance = a, l.farDistance = s), l
	}

	function E(i, r, a, s, o) {
		if (!1 === i.visible) return;
		if (i.layers.test(r.layers) && (i.isMesh || i.isLine || i.isPoints) && (i.castShadow || i.receiveShadow && 3 === o) && (!i.frustumCulled || n.intersectsObject(i))) {
			i.modelViewMatrix.multiplyMatrices(a.matrixWorldInverse, i.matrixWorld);
			const n = t.update(i),
				r = i.material;
			if (Array.isArray(r)) {
				const t = n.groups;
				for (let l = 0, c = t.length; l < c; l++) {
					const c = t[l],
						h = r[c.materialIndex];
					if (h && h.visible) {
						const t = y(i, n, h, s, a.near, a.far, o);
						e.renderBufferDirect(a, null, n, t, i, c)
					}
				}
			} else if (r.visible) {
				const t = y(i, n, r, s, a.near, a.far, o);
				e.renderBufferDirect(a, null, n, t, i, null)
			}
		}
		const l = i.children;
		for (let e = 0, t = l.length; e < t; e++) E(l[e], r, a, s, o)
	}
	this.enabled = !1, this.autoUpdate = !0, this.needsUpdate = !1, this.type = 1, this.render = function(t, o, l) {
		if (!1 === A.enabled) return;
		if (!1 === A.autoUpdate && !1 === A.needsUpdate) return;
		if (0 === t.length) return;
		const c = e.getRenderTarget(),
			h = e.getActiveCubeFace(),
			u = e.getActiveMipmapLevel(),
			d = e.state;
		d.setBlending(0), d.buffers.color.setClear(1, 1, 1, 1), d.buffers.depth.setTest(!0), d.setScissorTest(!1);
		for (let c = 0, h = t.length; c < h; c++) {
			const h = t[c],
				u = h.shadow;
			if (void 0 === u) {
				console.warn("THREE.WebGLShadowMap:", h, "has no shadow.");
				continue
			}
			if (!1 === u.autoUpdate && !1 === u.needsUpdate) continue;
			r.copy(u.mapSize);
			const p = u.getFrameExtents();
			if (r.multiply(p), a.copy(u.mapSize), (r.x > i || r.y > i) && (r.x > i && (a.x = Math.floor(i / p.x), r.x = a.x * p.x, u.mapSize.x = a.x), r.y > i && (a.y = Math.floor(i / p.y), r.y = a.y * p.y, u.mapSize.y = a.y)), null === u.map && !u.isPointLightShadow && 3 === this.type) {
				const e = {
					minFilter: 1006,
					magFilter: 1006,
					format: 1023
				};
				u.map = new WebGLRenderTarget(r.x, r.y, e), u.map.texture.name = h.name + ".shadowMap", u.mapPass = new WebGLRenderTarget(r.x, r.y, e), u.camera.updateProjectionMatrix()
			}
			if (null === u.map) {
				const e = {
					minFilter: 1003,
					magFilter: 1003,
					format: 1023
				};
				u.map = new WebGLRenderTarget(r.x, r.y, e), u.map.texture.name = h.name + ".shadowMap", u.camera.updateProjectionMatrix()
			}
			e.setRenderTarget(u.map), e.clear();
			const m = u.getViewportCount();
			for (let e = 0; e < m; e++) {
				const t = u.getViewport(e);
				s.set(a.x * t.x, a.y * t.y, a.x * t.z, a.y * t.w), d.viewport(s), u.updateMatrices(h, e), n = u.getFrustum(), E(o, l, u.camera, h, this.type)
			}
			u.isPointLightShadow || 3 !== this.type || g(u, l), u.needsUpdate = !1
		}
		A.needsUpdate = !1, e.setRenderTarget(c, h, u)
	}
}

function WebGLState(e, t, i) {
	const n = i.isWebGL2;
	const r = new function() {
			let t = !1;
			const i = new Vector4;
			let n = null;
			const r = new Vector4(0, 0, 0, 0);
			return {
				setMask: function(i) {
					n === i || t || (e.colorMask(i, i, i, i), n = i)
				},
				setLocked: function(e) {
					t = e
				},
				setClear: function(t, n, a, s, o) {
					!0 === o && (t *= s, n *= s, a *= s), i.set(t, n, a, s), !1 === r.equals(i) && (e.clearColor(t, n, a, s), r.copy(i))
				},
				reset: function() {
					t = !1, n = null, r.set(-1, 0, 0, 0)
				}
			}
		},
		a = new function() {
			let t = !1,
				i = null,
				n = null,
				r = null;
			return {
				setTest: function(e) {
					e ? D(2929) : P(2929)
				},
				setMask: function(n) {
					i === n || t || (e.depthMask(n), i = n)
				},
				setFunc: function(t) {
					if (n !== t) {
						if (t) switch (t) {
							case 0:
								e.depthFunc(512);
								break;
							case 1:
								e.depthFunc(519);
								break;
							case 2:
								e.depthFunc(513);
								break;
							case 3:
								e.depthFunc(515);
								break;
							case 4:
								e.depthFunc(514);
								break;
							case 5:
								e.depthFunc(518);
								break;
							case 6:
								e.depthFunc(516);
								break;
							case 7:
								e.depthFunc(517);
								break;
							default:
								e.depthFunc(515)
						} else e.depthFunc(515);
						n = t
					}
				},
				setLocked: function(e) {
					t = e
				},
				setClear: function(t) {
					r !== t && (e.clearDepth(t), r = t)
				},
				reset: function() {
					t = !1, i = null, n = null, r = null
				}
			}
		},
		s = new function() {
			let t = !1,
				i = null,
				n = null,
				r = null,
				a = null,
				s = null,
				o = null,
				l = null,
				c = null;
			return {
				setTest: function(e) {
					t || (e ? D(2960) : P(2960))
				},
				setMask: function(n) {
					i === n || t || (e.stencilMask(n), i = n)
				},
				setFunc: function(t, i, s) {
					n === t && r === i && a === s || (e.stencilFunc(t, i, s), n = t, r = i, a = s)
				},
				setOp: function(t, i, n) {
					s === t && o === i && l === n || (e.stencilOp(t, i, n), s = t, o = i, l = n)
				},
				setLocked: function(e) {
					t = e
				},
				setClear: function(t) {
					c !== t && (e.clearStencil(t), c = t)
				},
				reset: function() {
					t = !1, i = null, n = null, r = null, a = null, s = null, o = null, l = null, c = null
				}
			}
		};
	let o = {},
		l = null,
		c = !1,
		h = null,
		u = null,
		d = null,
		p = null,
		m = null,
		A = null,
		g = null,
		f = !1,
		v = null,
		y = null,
		E = null,
		_ = null,
		b = null;
	const x = e.getParameter(35661);
	let w = !1,
		C = 0;
	const S = e.getParameter(7938); - 1 !== S.indexOf("WebGL") ? (C = parseFloat(/^WebGL (\d)/.exec(S)[1]), w = C >= 1) : -1 !== S.indexOf("OpenGL ES") && (C = parseFloat(/^OpenGL ES (\d)/.exec(S)[1]), w = C >= 2);
	let I = null,
		M = {};
	const T = new Vector4,
		B = new Vector4;

	function L(t, i, n) {
		const r = new Uint8Array(4),
			a = e.createTexture();
		e.bindTexture(t, a), e.texParameteri(t, 10241, 9728), e.texParameteri(t, 10240, 9728);
		for (let t = 0; t < n; t++) e.texImage2D(i + t, 0, 6408, 1, 1, 0, 6408, 5121, r);
		return a
	}
	const R = {};

	function D(t) {
		!0 !== o[t] && (e.enable(t), o[t] = !0)
	}

	function P(t) {
		!1 !== o[t] && (e.disable(t), o[t] = !1)
	}
	R[3553] = L(3553, 3553, 1), R[34067] = L(34067, 34069, 6), r.setClear(0, 0, 0, 1), a.setClear(1), s.setClear(0), D(2929), a.setFunc(3), N(!1), k(1), D(2884), O(0);
	const Q = {
		100: 32774,
		101: 32778,
		102: 32779
	};
	if (n) Q[103] = 32775, Q[104] = 32776;
	else {
		const e = t.get("EXT_blend_minmax");
		null !== e && (Q[103] = e.MIN_EXT, Q[104] = e.MAX_EXT)
	}
	const F = {
		200: 0,
		201: 1,
		202: 768,
		204: 770,
		210: 776,
		208: 774,
		206: 772,
		203: 769,
		205: 771,
		209: 775,
		207: 773
	};

	function O(t, i, n, r, a, s, o, l) {
		if (0 !== t) {
			if (!1 === c && (D(3042), c = !0), 5 === t) a = a || i, s = s || n, o = o || r, i === u && a === m || (e.blendEquationSeparate(Q[i], Q[a]), u = i, m = a), n === d && r === p && s === A && o === g || (e.blendFuncSeparate(F[n], F[r], F[s], F[o]), d = n, p = r, A = s, g = o), h = t, f = null;
			else if (t !== h || l !== f) {
				if (100 === u && 100 === m || (e.blendEquation(32774), u = 100, m = 100), l) switch (t) {
					case 1:
						e.blendFuncSeparate(1, 771, 1, 771);
						break;
					case 2:
						e.blendFunc(1, 1);
						break;
					case 3:
						e.blendFuncSeparate(0, 0, 769, 771);
						break;
					case 4:
						e.blendFuncSeparate(0, 768, 0, 770);
						break;
					default:
						console.error("THREE.WebGLState: Invalid blending: ", t)
				} else switch (t) {
					case 1:
						e.blendFuncSeparate(770, 771, 1, 771);
						break;
					case 2:
						e.blendFunc(770, 1);
						break;
					case 3:
						e.blendFunc(0, 769);
						break;
					case 4:
						e.blendFunc(0, 768);
						break;
					default:
						console.error("THREE.WebGLState: Invalid blending: ", t)
				}
				d = null, p = null, A = null, g = null, h = t, f = l
			}
		} else !0 === c && (P(3042), c = !1)
	}

	function N(t) {
		v !== t && (t ? e.frontFace(2304) : e.frontFace(2305), v = t)
	}

	function k(t) {
		0 !== t ? (D(2884), t !== y && (1 === t ? e.cullFace(1029) : 2 === t ? e.cullFace(1028) : e.cullFace(1032))) : P(2884), y = t
	}

	function U(t, i, n) {
		t ? (D(32823), _ === i && b === n || (e.polygonOffset(i, n), _ = i, b = n)) : P(32823)
	}

	function G(t) {
		void 0 === t && (t = 33984 + x - 1), I !== t && (e.activeTexture(t), I = t)
	}
	return {
		buffers: {
			color: r,
			depth: a,
			stencil: s
		},
		enable: D,
		disable: P,
		useProgram: function(t) {
			return l !== t && (e.useProgram(t), l = t, !0)
		},
		setBlending: O,
		setMaterial: function(e, t) {
			2 === e.side ? P(2884) : D(2884);
			let i = 1 === e.side;
			t && (i = !i), N(i), 1 === e.blending && !1 === e.transparent ? O(0) : O(e.blending, e.blendEquation, e.blendSrc, e.blendDst, e.blendEquationAlpha, e.blendSrcAlpha, e.blendDstAlpha, e.premultipliedAlpha), a.setFunc(e.depthFunc), a.setTest(e.depthTest), a.setMask(e.depthWrite), r.setMask(e.colorWrite);
			const n = e.stencilWrite;
			s.setTest(n), n && (s.setMask(e.stencilWriteMask), s.setFunc(e.stencilFunc, e.stencilRef, e.stencilFuncMask), s.setOp(e.stencilFail, e.stencilZFail, e.stencilZPass)), U(e.polygonOffset, e.polygonOffsetFactor, e.polygonOffsetUnits)
		},
		setFlipSided: N,
		setCullFace: k,
		setLineWidth: function(t) {
			t !== E && (w && e.lineWidth(t), E = t)
		},
		setPolygonOffset: U,
		setScissorTest: function(e) {
			e ? D(3089) : P(3089)
		},
		activeTexture: G,
		bindTexture: function(t, i) {
			null === I && G();
			let n = M[I];
			void 0 === n && (n = {
				type: void 0,
				texture: void 0
			}, M[I] = n), n.type === t && n.texture === i || (e.bindTexture(t, i || R[t]), n.type = t, n.texture = i)
		},
		unbindTexture: function() {
			const t = M[I];
			void 0 !== t && void 0 !== t.type && (e.bindTexture(t.type, null), t.type = void 0, t.texture = void 0)
		},
		compressedTexImage2D: function() {
			try {
				e.compressedTexImage2D.apply(e, arguments)
			} catch (e) {
				console.error("THREE.WebGLState:", e)
			}
		},
		texImage2D: function() {
			try {
				e.texImage2D.apply(e, arguments)
			} catch (e) {
				console.error("THREE.WebGLState:", e)
			}
		},
		texImage3D: function() {
			try {
				e.texImage3D.apply(e, arguments)
			} catch (e) {
				console.error("THREE.WebGLState:", e)
			}
		},
		scissor: function(t) {
			!1 === T.equals(t) && (e.scissor(t.x, t.y, t.z, t.w), T.copy(t))
		},
		viewport: function(t) {
			!1 === B.equals(t) && (e.viewport(t.x, t.y, t.z, t.w), B.copy(t))
		},
		reset: function() {
			e.disable(3042), e.disable(2884), e.disable(2929), e.disable(32823), e.disable(3089), e.disable(2960), e.blendEquation(32774), e.blendFunc(1, 0), e.blendFuncSeparate(1, 0, 1, 0), e.colorMask(!0, !0, !0, !0), e.clearColor(0, 0, 0, 0), e.depthMask(!0), e.depthFunc(513), e.clearDepth(1), e.stencilMask(4294967295), e.stencilFunc(519, 0, 4294967295), e.stencilOp(7680, 7680, 7680), e.clearStencil(0), e.cullFace(1029), e.frontFace(2305), e.polygonOffset(0, 0), e.activeTexture(33984), e.useProgram(null), e.lineWidth(1), e.scissor(0, 0, e.canvas.width, e.canvas.height), e.viewport(0, 0, e.canvas.width, e.canvas.height), o = {}, I = null, M = {}, l = null, c = !1, h = null, u = null, d = null, p = null, m = null, A = null, g = null, f = !1, v = null, y = null, E = null, _ = null, b = null, r.reset(), a.reset(), s.reset()
		}
	}
}

function WebGLTextures(e, t, i, n, r, a, s) {
	const o = r.isWebGL2,
		l = r.maxTextures,
		c = r.maxCubemapSize,
		h = r.maxTextureSize,
		u = r.maxSamples,
		d = new WeakMap;
	let p, m = !1;
	try {
		m = "undefined" != typeof OffscreenCanvas && null !== new OffscreenCanvas(1, 1).getContext("2d")
	} catch (e) {}

	function A(e, t) {
		return m ? new OffscreenCanvas(e, t) : document.createElementNS("http://www.w3.org/1999/xhtml", "canvas")
	}

	function g(e, t, i, n) {
		let r = 1;
		if ((e.width > n || e.height > n) && (r = n / Math.max(e.width, e.height)), r < 1 || !0 === t) {
			if ("undefined" != typeof HTMLImageElement && e instanceof HTMLImageElement || "undefined" != typeof HTMLCanvasElement && e instanceof HTMLCanvasElement || "undefined" != typeof ImageBitmap && e instanceof ImageBitmap) {
				const n = t ? MathUtils.floorPowerOfTwo : Math.floor,
					a = n(r * e.width),
					s = n(r * e.height);
				void 0 === p && (p = A(a, s));
				const o = i ? A(a, s) : p;
				o.width = a, o.height = s;
				return o.getContext("2d").drawImage(e, 0, 0, a, s), console.warn("THREE.WebGLRenderer: Texture has been resized from (" + e.width + "x" + e.height + ") to (" + a + "x" + s + ")."), o
			}
			return "data" in e && console.warn("THREE.WebGLRenderer: Image in DataTexture is too big (" + e.width + "x" + e.height + ")."), e
		}
		return e
	}

	function f(e) {
		return MathUtils.isPowerOfTwo(e.width) && MathUtils.isPowerOfTwo(e.height)
	}

	function v(e, t) {
		return e.generateMipmaps && t && 1003 !== e.minFilter && 1006 !== e.minFilter
	}

	function y(t, i, r, a) {
		e.generateMipmap(t);
		n.get(i).__maxMipLevel = Math.log2(Math.max(r, a))
	}

	function E(i, n, r) {
		if (!1 === o) return n;
		if (null !== i) {
			if (void 0 !== e[i]) return e[i];
			console.warn("THREE.WebGLRenderer: Attempt to use non-existing WebGL internal format '" + i + "'")
		}
		let a = n;
		return 6403 === n && (5126 === r && (a = 33326), 5131 === r && (a = 33325), 5121 === r && (a = 33321)), 6407 === n && (5126 === r && (a = 34837), 5131 === r && (a = 34843), 5121 === r && (a = 32849)), 6408 === n && (5126 === r && (a = 34836), 5131 === r && (a = 34842), 5121 === r && (a = 32856)), 33325 !== a && 33326 !== a && 34842 !== a && 34836 !== a || t.get("EXT_color_buffer_float"), a
	}

	function _(e) {
		return 1003 === e || 1004 === e || 1005 === e ? 9728 : 9729
	}

	function b(t) {
		const i = t.target;
		i.removeEventListener("dispose", b),
			function(t) {
				const i = n.get(t);
				if (void 0 === i.__webglInit) return;
				e.deleteTexture(i.__webglTexture), n.remove(t)
			}(i), i.isVideoTexture && d.delete(i), s.memory.textures--
	}

	function x(t) {
		const i = t.target;
		i.removeEventListener("dispose", x),
			function(t) {
				const i = t.texture,
					r = n.get(t),
					a = n.get(i);
				if (!t) return;
				void 0 !== a.__webglTexture && e.deleteTexture(a.__webglTexture);
				t.depthTexture && t.depthTexture.dispose();
				if (t.isWebGLCubeRenderTarget)
					for (let t = 0; t < 6; t++) e.deleteFramebuffer(r.__webglFramebuffer[t]), r.__webglDepthbuffer && e.deleteRenderbuffer(r.__webglDepthbuffer[t]);
				else e.deleteFramebuffer(r.__webglFramebuffer), r.__webglDepthbuffer && e.deleteRenderbuffer(r.__webglDepthbuffer), r.__webglMultisampledFramebuffer && e.deleteFramebuffer(r.__webglMultisampledFramebuffer), r.__webglColorRenderbuffer && e.deleteRenderbuffer(r.__webglColorRenderbuffer), r.__webglDepthRenderbuffer && e.deleteRenderbuffer(r.__webglDepthRenderbuffer);
				n.remove(i), n.remove(t)
			}(i), s.memory.textures--
	}
	let w = 0;

	function C(e, t) {
		const r = n.get(e);
		if (e.isVideoTexture && function(e) {
				const t = s.render.frame;
				d.get(e) !== t && (d.set(e, t), e.update())
			}(e), e.version > 0 && r.__version !== e.version) {
			const i = e.image;
			if (void 0 === i) console.warn("THREE.WebGLRenderer: Texture marked for update but image is undefined");
			else {
				if (!1 !== i.complete) return void L(r, e, t);
				console.warn("THREE.WebGLRenderer: Texture marked for update but image is incomplete")
			}
		}
		i.activeTexture(33984 + t), i.bindTexture(3553, r.__webglTexture)
	}

	function S(t, r) {
		const s = n.get(t);
		t.version > 0 && s.__version !== t.version ? function(t, n, r) {
			if (6 !== n.image.length) return;
			B(t, n), i.activeTexture(33984 + r), i.bindTexture(34067, t.__webglTexture), e.pixelStorei(37440, n.flipY), e.pixelStorei(37441, n.premultiplyAlpha), e.pixelStorei(3317, n.unpackAlignment), e.pixelStorei(37443, 0);
			const s = n && (n.isCompressedTexture || n.image[0].isCompressedTexture),
				l = n.image[0] && n.image[0].isDataTexture,
				h = [];
			for (let e = 0; e < 6; e++) h[e] = s || l ? l ? n.image[e].image : n.image[e] : g(n.image[e], !1, !0, c);
			const u = h[0],
				d = f(u) || o,
				p = a.convert(n.format),
				m = a.convert(n.type),
				A = E(n.internalFormat, p, m);
			let _;
			if (T(34067, n, d), s) {
				for (let e = 0; e < 6; e++) {
					_ = h[e].mipmaps;
					for (let t = 0; t < _.length; t++) {
						const r = _[t];
						1023 !== n.format && 1022 !== n.format ? null !== p ? i.compressedTexImage2D(34069 + e, t, A, r.width, r.height, 0, r.data) : console.warn("THREE.WebGLRenderer: Attempt to load unsupported compressed texture format in .setTextureCube()") : i.texImage2D(34069 + e, t, A, r.width, r.height, 0, p, m, r.data)
					}
				}
				t.__maxMipLevel = _.length - 1
			} else {
				_ = n.mipmaps;
				for (let e = 0; e < 6; e++)
					if (l) {
						i.texImage2D(34069 + e, 0, A, h[e].width, h[e].height, 0, p, m, h[e].data);
						for (let t = 0; t < _.length; t++) {
							const n = _[t].image[e].image;
							i.texImage2D(34069 + e, t + 1, A, n.width, n.height, 0, p, m, n.data)
						}
					} else {
						i.texImage2D(34069 + e, 0, A, p, m, h[e]);
						for (let t = 0; t < _.length; t++) {
							const n = _[t];
							i.texImage2D(34069 + e, t + 1, A, p, m, n.image[e])
						}
					} t.__maxMipLevel = _.length
			}
			v(n, d) && y(34067, n, u.width, u.height);
			t.__version = n.version, n.onUpdate && n.onUpdate(n)
		}(s, t, r) : (i.activeTexture(33984 + r), i.bindTexture(34067, s.__webglTexture))
	}
	const I = {
			1e3: 10497,
			1001: 33071,
			1002: 33648
		},
		M = {
			1003: 9728,
			1004: 9984,
			1005: 9986,
			1006: 9729,
			1007: 9985,
			1008: 9987
		};

	function T(i, a, s) {
		if (s ? (e.texParameteri(i, 10242, I[a.wrapS]), e.texParameteri(i, 10243, I[a.wrapT]), 32879 !== i && 35866 !== i || e.texParameteri(i, 32882, I[a.wrapR]), e.texParameteri(i, 10240, M[a.magFilter]), e.texParameteri(i, 10241, M[a.minFilter])) : (e.texParameteri(i, 10242, 33071), e.texParameteri(i, 10243, 33071), 32879 !== i && 35866 !== i || e.texParameteri(i, 32882, 33071), 1001 === a.wrapS && 1001 === a.wrapT || console.warn("THREE.WebGLRenderer: Texture is not power of two. Texture.wrapS and Texture.wrapT should be set to THREE.ClampToEdgeWrapping."), e.texParameteri(i, 10240, _(a.magFilter)), e.texParameteri(i, 10241, _(a.minFilter)), 1003 !== a.minFilter && 1006 !== a.minFilter && console.warn("THREE.WebGLRenderer: Texture is not power of two. Texture.minFilter should be set to THREE.NearestFilter or THREE.LinearFilter.")), !0 === t.has("EXT_texture_filter_anisotropic")) {
			const s = t.get("EXT_texture_filter_anisotropic");
			if (1015 === a.type && !1 === t.has("OES_texture_float_linear")) return;
			if (!1 === o && 1016 === a.type && !1 === t.has("OES_texture_half_float_linear")) return;
			(a.anisotropy > 1 || n.get(a).__currentAnisotropy) && (e.texParameterf(i, s.TEXTURE_MAX_ANISOTROPY_EXT, Math.min(a.anisotropy, r.getMaxAnisotropy())), n.get(a).__currentAnisotropy = a.anisotropy)
		}
	}

	function B(t, i) {
		void 0 === t.__webglInit && (t.__webglInit = !0, i.addEventListener("dispose", b), t.__webglTexture = e.createTexture(), s.memory.textures++)
	}

	function L(t, n, r) {
		let s = 3553;
		n.isDataTexture2DArray && (s = 35866), n.isDataTexture3D && (s = 32879), B(t, n), i.activeTexture(33984 + r), i.bindTexture(s, t.__webglTexture), e.pixelStorei(37440, n.flipY), e.pixelStorei(37441, n.premultiplyAlpha), e.pixelStorei(3317, n.unpackAlignment), e.pixelStorei(37443, 0);
		const l = function(e) {
				return !o && (1001 !== e.wrapS || 1001 !== e.wrapT || 1003 !== e.minFilter && 1006 !== e.minFilter)
			}(n) && !1 === f(n.image),
			c = g(n.image, l, !1, h),
			u = f(c) || o,
			d = a.convert(n.format);
		let p, m = a.convert(n.type),
			A = E(n.internalFormat, d, m);
		T(s, n, u);
		const _ = n.mipmaps;
		if (n.isDepthTexture) A = 6402, o ? A = 1015 === n.type ? 36012 : 1014 === n.type ? 33190 : 1020 === n.type ? 35056 : 33189 : 1015 === n.type && console.error("WebGLRenderer: Floating point depth texture requires WebGL2."), 1026 === n.format && 6402 === A && 1012 !== n.type && 1014 !== n.type && (console.warn("THREE.WebGLRenderer: Use UnsignedShortType or UnsignedIntType for DepthFormat DepthTexture."), n.type = 1012, m = a.convert(n.type)), 1027 === n.format && 6402 === A && (A = 34041, 1020 !== n.type && (console.warn("THREE.WebGLRenderer: Use UnsignedInt248Type for DepthStencilFormat DepthTexture."), n.type = 1020, m = a.convert(n.type))), i.texImage2D(3553, 0, A, c.width, c.height, 0, d, m, null);
		else if (n.isDataTexture)
			if (_.length > 0 && u) {
				for (let e = 0, t = _.length; e < t; e++) p = _[e], i.texImage2D(3553, e, A, p.width, p.height, 0, d, m, p.data);
				n.generateMipmaps = !1, t.__maxMipLevel = _.length - 1
			} else i.texImage2D(3553, 0, A, c.width, c.height, 0, d, m, c.data), t.__maxMipLevel = 0;
		else if (n.isCompressedTexture) {
			for (let e = 0, t = _.length; e < t; e++) p = _[e], 1023 !== n.format && 1022 !== n.format ? null !== d ? i.compressedTexImage2D(3553, e, A, p.width, p.height, 0, p.data) : console.warn("THREE.WebGLRenderer: Attempt to load unsupported compressed texture format in .uploadTexture()") : i.texImage2D(3553, e, A, p.width, p.height, 0, d, m, p.data);
			t.__maxMipLevel = _.length - 1
		} else if (n.isDataTexture2DArray) i.texImage3D(35866, 0, A, c.width, c.height, c.depth, 0, d, m, c.data), t.__maxMipLevel = 0;
		else if (n.isDataTexture3D) i.texImage3D(32879, 0, A, c.width, c.height, c.depth, 0, d, m, c.data), t.__maxMipLevel = 0;
		else if (_.length > 0 && u) {
			for (let e = 0, t = _.length; e < t; e++) p = _[e], i.texImage2D(3553, e, A, d, m, p);
			n.generateMipmaps = !1, t.__maxMipLevel = _.length - 1
		} else i.texImage2D(3553, 0, A, d, m, c), t.__maxMipLevel = 0;
		v(n, u) && y(s, n, c.width, c.height), t.__version = n.version, n.onUpdate && n.onUpdate(n)
	}

	function R(t, r, s, o) {
		const l = r.texture,
			c = a.convert(l.format),
			h = a.convert(l.type),
			u = E(l.internalFormat, c, h);
		32879 === o || 35866 === o ? i.texImage3D(o, 0, u, r.width, r.height, r.depth, 0, c, h, null) : i.texImage2D(o, 0, u, r.width, r.height, 0, c, h, null), e.bindFramebuffer(36160, t), e.framebufferTexture2D(36160, s, o, n.get(l).__webglTexture, 0), e.bindFramebuffer(36160, null)
	}

	function D(t, i, n) {
		if (e.bindRenderbuffer(36161, t), i.depthBuffer && !i.stencilBuffer) {
			let r = 33189;
			if (n) {
				const t = i.depthTexture;
				t && t.isDepthTexture && (1015 === t.type ? r = 36012 : 1014 === t.type && (r = 33190));
				const n = Q(i);
				e.renderbufferStorageMultisample(36161, n, r, i.width, i.height)
			} else e.renderbufferStorage(36161, r, i.width, i.height);
			e.framebufferRenderbuffer(36160, 36096, 36161, t)
		} else if (i.depthBuffer && i.stencilBuffer) {
			if (n) {
				const t = Q(i);
				e.renderbufferStorageMultisample(36161, t, 35056, i.width, i.height)
			} else e.renderbufferStorage(36161, 34041, i.width, i.height);
			e.framebufferRenderbuffer(36160, 33306, 36161, t)
		} else {
			const t = i.texture,
				r = a.convert(t.format),
				s = a.convert(t.type),
				o = E(t.internalFormat, r, s);
			if (n) {
				const t = Q(i);
				e.renderbufferStorageMultisample(36161, t, o, i.width, i.height)
			} else e.renderbufferStorage(36161, o, i.width, i.height)
		}
		e.bindRenderbuffer(36161, null)
	}

	function P(t) {
		const i = n.get(t),
			r = !0 === t.isWebGLCubeRenderTarget;
		if (t.depthTexture) {
			if (r) throw new Error("target.depthTexture not supported in Cube render targets");
			! function(t, i) {
				if (i && i.isWebGLCubeRenderTarget) throw new Error("Depth Texture with cube render targets is not supported");
				if (e.bindFramebuffer(36160, t), !i.depthTexture || !i.depthTexture.isDepthTexture) throw new Error("renderTarget.depthTexture must be an instance of THREE.DepthTexture");
				n.get(i.depthTexture).__webglTexture && i.depthTexture.image.width === i.width && i.depthTexture.image.height === i.height || (i.depthTexture.image.width = i.width, i.depthTexture.image.height = i.height, i.depthTexture.needsUpdate = !0), C(i.depthTexture, 0);
				const r = n.get(i.depthTexture).__webglTexture;
				if (1026 === i.depthTexture.format) e.framebufferTexture2D(36160, 36096, 3553, r, 0);
				else {
					if (1027 !== i.depthTexture.format) throw new Error("Unknown depthTexture format");
					e.framebufferTexture2D(36160, 33306, 3553, r, 0)
				}
			}(i.__webglFramebuffer, t)
		} else if (r) {
			i.__webglDepthbuffer = [];
			for (let n = 0; n < 6; n++) e.bindFramebuffer(36160, i.__webglFramebuffer[n]), i.__webglDepthbuffer[n] = e.createRenderbuffer(), D(i.__webglDepthbuffer[n], t, !1)
		} else e.bindFramebuffer(36160, i.__webglFramebuffer), i.__webglDepthbuffer = e.createRenderbuffer(), D(i.__webglDepthbuffer, t, !1);
		e.bindFramebuffer(36160, null)
	}

	function Q(e) {
		return o && e.isWebGLMultisampleRenderTarget ? Math.min(u, e.samples) : 0
	}
	let F = !1,
		O = !1;
	this.allocateTextureUnit = function() {
		const e = w;
		return e >= l && console.warn("THREE.WebGLTextures: Trying to use " + e + " texture units while this GPU supports only " + l), w += 1, e
	}, this.resetTextureUnits = function() {
		w = 0
	}, this.setTexture2D = C, this.setTexture2DArray = function(e, t) {
		const r = n.get(e);
		e.version > 0 && r.__version !== e.version ? L(r, e, t) : (i.activeTexture(33984 + t), i.bindTexture(35866, r.__webglTexture))
	}, this.setTexture3D = function(e, t) {
		const r = n.get(e);
		e.version > 0 && r.__version !== e.version ? L(r, e, t) : (i.activeTexture(33984 + t), i.bindTexture(32879, r.__webglTexture))
	}, this.setTextureCube = S, this.setupRenderTarget = function(t) {
		const r = t.texture,
			l = n.get(t),
			c = n.get(r);
		t.addEventListener("dispose", x), c.__webglTexture = e.createTexture(), s.memory.textures++;
		const h = !0 === t.isWebGLCubeRenderTarget,
			u = !0 === t.isWebGLMultisampleRenderTarget,
			d = r.isDataTexture3D || r.isDataTexture2DArray,
			p = f(t) || o;
		if (!o || 1022 !== r.format || 1015 !== r.type && 1016 !== r.type || (r.format = 1023, console.warn("THREE.WebGLRenderer: Rendering to textures with RGB format is not supported. Using RGBA format instead.")), h) {
			l.__webglFramebuffer = [];
			for (let t = 0; t < 6; t++) l.__webglFramebuffer[t] = e.createFramebuffer()
		} else if (l.__webglFramebuffer = e.createFramebuffer(), u)
			if (o) {
				l.__webglMultisampledFramebuffer = e.createFramebuffer(), l.__webglColorRenderbuffer = e.createRenderbuffer(), e.bindRenderbuffer(36161, l.__webglColorRenderbuffer);
				const i = a.convert(r.format),
					n = a.convert(r.type),
					s = E(r.internalFormat, i, n),
					o = Q(t);
				e.renderbufferStorageMultisample(36161, o, s, t.width, t.height), e.bindFramebuffer(36160, l.__webglMultisampledFramebuffer), e.framebufferRenderbuffer(36160, 36064, 36161, l.__webglColorRenderbuffer), e.bindRenderbuffer(36161, null), t.depthBuffer && (l.__webglDepthRenderbuffer = e.createRenderbuffer(), D(l.__webglDepthRenderbuffer, t, !0)), e.bindFramebuffer(36160, null)
			} else console.warn("THREE.WebGLRenderer: WebGLMultisampleRenderTarget can only be used with WebGL2.");
		if (h) {
			i.bindTexture(34067, c.__webglTexture), T(34067, r, p);
			for (let e = 0; e < 6; e++) R(l.__webglFramebuffer[e], t, 36064, 34069 + e);
			v(r, p) && y(34067, r, t.width, t.height), i.bindTexture(34067, null)
		} else {
			let e = 3553;
			if (d)
				if (o) {
					e = r.isDataTexture3D ? 32879 : 35866
				} else console.warn("THREE.DataTexture3D and THREE.DataTexture2DArray only supported with WebGL2.");
			i.bindTexture(e, c.__webglTexture), T(e, r, p), R(l.__webglFramebuffer, t, 36064, e), v(r, p) && y(3553, r, t.width, t.height), i.bindTexture(3553, null)
		}
		t.depthBuffer && P(t)
	}, this.updateRenderTargetMipmap = function(e) {
		const t = e.texture;
		if (v(t, f(e) || o)) {
			const r = e.isWebGLCubeRenderTarget ? 34067 : 3553,
				a = n.get(t).__webglTexture;
			i.bindTexture(r, a), y(r, t, e.width, e.height), i.bindTexture(r, null)
		}
	}, this.updateMultisampleRenderTarget = function(t) {
		if (t.isWebGLMultisampleRenderTarget)
			if (o) {
				const i = n.get(t);
				e.bindFramebuffer(36008, i.__webglMultisampledFramebuffer), e.bindFramebuffer(36009, i.__webglFramebuffer);
				const r = t.width,
					a = t.height;
				let s = 16384;
				t.depthBuffer && (s |= 256), t.stencilBuffer && (s |= 1024), e.blitFramebuffer(0, 0, r, a, 0, 0, r, a, s, 9728), e.bindFramebuffer(36160, i.__webglMultisampledFramebuffer)
			} else console.warn("THREE.WebGLRenderer: WebGLMultisampleRenderTarget can only be used with WebGL2.")
	}, this.safeSetTexture2D = function(e, t) {
		e && e.isWebGLRenderTarget && (!1 === F && (console.warn("THREE.WebGLTextures.safeSetTexture2D: don't use render targets as textures. Use their .texture property instead."), F = !0), e = e.texture), C(e, t)
	}, this.safeSetTextureCube = function(e, t) {
		e && e.isWebGLCubeRenderTarget && (!1 === O && (console.warn("THREE.WebGLTextures.safeSetTextureCube: don't use cube render targets as textures. Use their .texture property instead."), O = !0), e = e.texture), S(e, t)
	}
}

function WebGLUtils(e, t, i) {
	const n = i.isWebGL2;
	return {
		convert: function(e) {
			let i;
			if (1009 === e) return 5121;
			if (1017 === e) return 32819;
			if (1018 === e) return 32820;
			if (1019 === e) return 33635;
			if (1010 === e) return 5120;
			if (1011 === e) return 5122;
			if (1012 === e) return 5123;
			if (1013 === e) return 5124;
			if (1014 === e) return 5125;
			if (1015 === e) return 5126;
			if (1016 === e) return n ? 5131 : (i = t.get("OES_texture_half_float"), null !== i ? i.HALF_FLOAT_OES : null);
			if (1021 === e) return 6406;
			if (1022 === e) return 6407;
			if (1023 === e) return 6408;
			if (1024 === e) return 6409;
			if (1025 === e) return 6410;
			if (1026 === e) return 6402;
			if (1027 === e) return 34041;
			if (1028 === e) return 6403;
			if (1029 === e) return 36244;
			if (1030 === e) return 33319;
			if (1031 === e) return 33320;
			if (1032 === e) return 36248;
			if (1033 === e) return 36249;
			if (33776 === e || 33777 === e || 33778 === e || 33779 === e) {
				if (i = t.get("WEBGL_compressed_texture_s3tc"), null === i) return null;
				if (33776 === e) return i.COMPRESSED_RGB_S3TC_DXT1_EXT;
				if (33777 === e) return i.COMPRESSED_RGBA_S3TC_DXT1_EXT;
				if (33778 === e) return i.COMPRESSED_RGBA_S3TC_DXT3_EXT;
				if (33779 === e) return i.COMPRESSED_RGBA_S3TC_DXT5_EXT
			}
			if (35840 === e || 35841 === e || 35842 === e || 35843 === e) {
				if (i = t.get("WEBGL_compressed_texture_pvrtc"), null === i) return null;
				if (35840 === e) return i.COMPRESSED_RGB_PVRTC_4BPPV1_IMG;
				if (35841 === e) return i.COMPRESSED_RGB_PVRTC_2BPPV1_IMG;
				if (35842 === e) return i.COMPRESSED_RGBA_PVRTC_4BPPV1_IMG;
				if (35843 === e) return i.COMPRESSED_RGBA_PVRTC_2BPPV1_IMG
			}
			if (36196 === e) return i = t.get("WEBGL_compressed_texture_etc1"), null !== i ? i.COMPRESSED_RGB_ETC1_WEBGL : null;
			if ((37492 === e || 37496 === e) && (i = t.get("WEBGL_compressed_texture_etc"), null !== i)) {
				if (37492 === e) return i.COMPRESSED_RGB8_ETC2;
				if (37496 === e) return i.COMPRESSED_RGBA8_ETC2_EAC
			}
			return 37808 === e || 37809 === e || 37810 === e || 37811 === e || 37812 === e || 37813 === e || 37814 === e || 37815 === e || 37816 === e || 37817 === e || 37818 === e || 37819 === e || 37820 === e || 37821 === e || 37840 === e || 37841 === e || 37842 === e || 37843 === e || 37844 === e || 37845 === e || 37846 === e || 37847 === e || 37848 === e || 37849 === e || 37850 === e || 37851 === e || 37852 === e || 37853 === e ? (i = t.get("WEBGL_compressed_texture_astc"), null !== i ? e : null) : 36492 === e ? (i = t.get("EXT_texture_compression_bptc"), null !== i ? e : null) : 1020 === e ? n ? 34042 : (i = t.get("WEBGL_depth_texture"), null !== i ? i.UNSIGNED_INT_24_8_WEBGL : null) : void 0
		}
	}
}

function ArrayCamera(e = []) {
	PerspectiveCamera.call(this), this.cameras = e
}
ArrayCamera.prototype = Object.assign(Object.create(PerspectiveCamera.prototype), {
	constructor: ArrayCamera,
	isArrayCamera: !0
});
class Group extends Object3D {
	constructor() {
		super(), this.type = "Group"
	}
}

function WebXRController() {
	this._targetRay = null, this._grip = null, this._hand = null
}

function WebXRManager(e, t) {
	const i = this;
	let n = null,
		r = 1,
		a = null,
		s = "local-floor",
		o = null;
	const l = [],
		c = new Map,
		h = new PerspectiveCamera;
	h.layers.enable(1), h.viewport = new Vector4;
	const u = new PerspectiveCamera;
	u.layers.enable(2), u.viewport = new Vector4;
	const d = [h, u],
		p = new ArrayCamera;
	p.layers.enable(1), p.layers.enable(2);
	let m = null,
		A = null;

	function g(e) {
		const t = c.get(e.inputSource);
		t && t.dispatchEvent({
			type: e.type,
			data: e.inputSource
		})
	}

	function f() {
		c.forEach((function(e, t) {
			e.disconnect(t)
		})), c.clear(), m = null, A = null, e.setFramebuffer(null), e.setRenderTarget(e.getRenderTarget()), x.stop(), i.isPresenting = !1, i.dispatchEvent({
			type: "sessionend"
		})
	}

	function v(e) {
		const t = n.inputSources;
		for (let e = 0; e < l.length; e++) c.set(t[e], l[e]);
		for (let t = 0; t < e.removed.length; t++) {
			const i = e.removed[t],
				n = c.get(i);
			n && (n.dispatchEvent({
				type: "disconnected",
				data: i
			}), c.delete(i))
		}
		for (let t = 0; t < e.added.length; t++) {
			const i = e.added[t],
				n = c.get(i);
			n && n.dispatchEvent({
				type: "connected",
				data: i
			})
		}
	}
	this.enabled = !1, this.isPresenting = !1, this.getController = function(e) {
		let t = l[e];
		return void 0 === t && (t = new WebXRController, l[e] = t), t.getTargetRaySpace()
	}, this.getControllerGrip = function(e) {
		let t = l[e];
		return void 0 === t && (t = new WebXRController, l[e] = t), t.getGripSpace()
	}, this.getHand = function(e) {
		let t = l[e];
		return void 0 === t && (t = new WebXRController, l[e] = t), t.getHandSpace()
	}, this.setFramebufferScaleFactor = function(e) {
		r = e, !0 === i.isPresenting && console.warn("THREE.WebXRManager: Cannot change framebuffer scale while presenting.")
	}, this.setReferenceSpaceType = function(e) {
		s = e, !0 === i.isPresenting && console.warn("THREE.WebXRManager: Cannot change reference space type while presenting.")
	}, this.getReferenceSpace = function() {
		return a
	}, this.getSession = function() {
		return n
	}, this.setSession = async function(e) {
		if (n = e, null !== n) {
			n.addEventListener("select", g), n.addEventListener("selectstart", g), n.addEventListener("selectend", g), n.addEventListener("squeeze", g), n.addEventListener("squeezestart", g), n.addEventListener("squeezeend", g), n.addEventListener("end", f), n.addEventListener("inputsourceschange", v);
			const e = t.getContextAttributes();
			!0 !== e.xrCompatible && await t.makeXRCompatible();
			const o = {
					antialias: e.antialias,
					alpha: e.alpha,
					depth: e.depth,
					stencil: e.stencil,
					framebufferScaleFactor: r
				},
				l = new XRWebGLLayer(n, t, o);
			n.updateRenderState({
				baseLayer: l
			}), a = await n.requestReferenceSpace(s), x.setContext(n), x.start(), i.isPresenting = !0, i.dispatchEvent({
				type: "sessionstart"
			})
		}
	};
	const y = new Vector3,
		E = new Vector3;

	function _(e, t) {
		null === t ? e.matrixWorld.copy(e.matrix) : e.matrixWorld.multiplyMatrices(t.matrixWorld, e.matrix), e.matrixWorldInverse.copy(e.matrixWorld).invert()
	}
	this.getCamera = function(e) {
		p.near = u.near = h.near = e.near, p.far = u.far = h.far = e.far, m === p.near && A === p.far || (n.updateRenderState({
			depthNear: p.near,
			depthFar: p.far
		}), m = p.near, A = p.far);
		const t = e.parent,
			i = p.cameras;
		_(p, t);
		for (let e = 0; e < i.length; e++) _(i[e], t);
		e.matrixWorld.copy(p.matrixWorld), e.matrix.copy(p.matrix), e.matrix.decompose(e.position, e.quaternion, e.scale);
		const r = e.children;
		for (let e = 0, t = r.length; e < t; e++) r[e].updateMatrixWorld(!0);
		return 2 === i.length ? function(e, t, i) {
			y.setFromMatrixPosition(t.matrixWorld), E.setFromMatrixPosition(i.matrixWorld);
			const n = y.distanceTo(E),
				r = t.projectionMatrix.elements,
				a = i.projectionMatrix.elements,
				s = r[14] / (r[10] - 1),
				o = r[14] / (r[10] + 1),
				l = (r[9] + 1) / r[5],
				c = (r[9] - 1) / r[5],
				h = (r[8] - 1) / r[0],
				u = (a[8] + 1) / a[0],
				d = s * h,
				p = s * u,
				m = n / (-h + u),
				A = m * -h;
			t.matrixWorld.decompose(e.position, e.quaternion, e.scale), e.translateX(A), e.translateZ(m), e.matrixWorld.compose(e.position, e.quaternion, e.scale), e.matrixWorldInverse.copy(e.matrixWorld).invert();
			const g = s + m,
				f = o + m,
				v = d - A,
				_ = p + (n - A),
				b = l * o / f * g,
				x = c * o / f * g;
			e.projectionMatrix.makePerspective(v, _, b, x, g, f)
		}(p, h, u) : p.projectionMatrix.copy(h.projectionMatrix), p
	};
	let b = null;
	const x = new WebGLAnimation;
	x.setAnimationLoop((function(t, i) {
		if (o = i.getViewerPose(a), null !== o) {
			const t = o.views,
				i = n.renderState.baseLayer;
			e.setFramebuffer(i.framebuffer);
			let r = !1;
			t.length !== p.cameras.length && (p.cameras.length = 0, r = !0);
			for (let e = 0; e < t.length; e++) {
				const n = t[e],
					a = i.getViewport(n),
					s = d[e];
				s.matrix.fromArray(n.transform.matrix), s.projectionMatrix.fromArray(n.projectionMatrix), s.viewport.set(a.x, a.y, a.width, a.height), 0 === e && p.matrix.copy(s.matrix), !0 === r && p.cameras.push(s)
			}
		}
		const r = n.inputSources;
		for (let e = 0; e < l.length; e++) {
			const t = l[e],
				n = r[e];
			t.update(n, i, a)
		}
		b && b(t, i)
	})), this.setAnimationLoop = function(e) {
		b = e
	}, this.dispose = function() {}
}

function WebGLMaterials(e) {
	function t(t, i) {
		t.opacity.value = i.opacity, i.color && t.diffuse.value.copy(i.color), i.emissive && t.emissive.value.copy(i.emissive).multiplyScalar(i.emissiveIntensity), i.map && (t.map.value = i.map), i.alphaMap && (t.alphaMap.value = i.alphaMap), i.specularMap && (t.specularMap.value = i.specularMap);
		const n = e.get(i).envMap;
		if (n) {
			t.envMap.value = n, t.flipEnvMap.value = n.isCubeTexture && n._needsFlipEnvMap ? -1 : 1, t.reflectivity.value = i.reflectivity, t.refractionRatio.value = i.refractionRatio;
			const r = e.get(n).__maxMipLevel;
			void 0 !== r && (t.maxMipLevel.value = r)
		}
		let r, a;
		i.lightMap && (t.lightMap.value = i.lightMap, t.lightMapIntensity.value = i.lightMapIntensity), i.aoMap && (t.aoMap.value = i.aoMap, t.aoMapIntensity.value = i.aoMapIntensity), i.map ? r = i.map : i.specularMap ? r = i.specularMap : i.displacementMap ? r = i.displacementMap : i.normalMap ? r = i.normalMap : i.bumpMap ? r = i.bumpMap : i.roughnessMap ? r = i.roughnessMap : i.metalnessMap ? r = i.metalnessMap : i.alphaMap ? r = i.alphaMap : i.emissiveMap ? r = i.emissiveMap : i.clearcoatMap ? r = i.clearcoatMap : i.clearcoatNormalMap ? r = i.clearcoatNormalMap : i.clearcoatRoughnessMap && (r = i.clearcoatRoughnessMap), void 0 !== r && (r.isWebGLRenderTarget && (r = r.texture), !0 === r.matrixAutoUpdate && r.updateMatrix(), t.uvTransform.value.copy(r.matrix)), i.aoMap ? a = i.aoMap : i.lightMap && (a = i.lightMap), void 0 !== a && (a.isWebGLRenderTarget && (a = a.texture), !0 === a.matrixAutoUpdate && a.updateMatrix(), t.uv2Transform.value.copy(a.matrix))
	}

	function i(t, i) {
		t.roughness.value = i.roughness, t.metalness.value = i.metalness, i.roughnessMap && (t.roughnessMap.value = i.roughnessMap), i.metalnessMap && (t.metalnessMap.value = i.metalnessMap), i.emissiveMap && (t.emissiveMap.value = i.emissiveMap), i.bumpMap && (t.bumpMap.value = i.bumpMap, t.bumpScale.value = i.bumpScale, 1 === i.side && (t.bumpScale.value *= -1)), i.normalMap && (t.normalMap.value = i.normalMap, t.normalScale.value.copy(i.normalScale), 1 === i.side && t.normalScale.value.negate()), i.displacementMap && (t.displacementMap.value = i.displacementMap, t.displacementScale.value = i.displacementScale, t.displacementBias.value = i.displacementBias);
		e.get(i).envMap && (t.envMapIntensity.value = i.envMapIntensity)
	}
	return {
		refreshFogUniforms: function(e, t) {
			e.fogColor.value.copy(t.color), t.isFog ? (e.fogNear.value = t.near, e.fogFar.value = t.far) : t.isFogExp2 && (e.fogDensity.value = t.density)
		},
		refreshMaterialUniforms: function(e, n, r, a) {
			n.isMeshBasicMaterial ? t(e, n) : n.isMeshLambertMaterial ? (t(e, n), function(e, t) {
				t.emissiveMap && (e.emissiveMap.value = t.emissiveMap)
			}(e, n)) : n.isMeshToonMaterial ? (t(e, n), function(e, t) {
				t.gradientMap && (e.gradientMap.value = t.gradientMap);
				t.emissiveMap && (e.emissiveMap.value = t.emissiveMap);
				t.bumpMap && (e.bumpMap.value = t.bumpMap, e.bumpScale.value = t.bumpScale, 1 === t.side && (e.bumpScale.value *= -1));
				t.normalMap && (e.normalMap.value = t.normalMap, e.normalScale.value.copy(t.normalScale), 1 === t.side && e.normalScale.value.negate());
				t.displacementMap && (e.displacementMap.value = t.displacementMap, e.displacementScale.value = t.displacementScale, e.displacementBias.value = t.displacementBias)
			}(e, n)) : n.isMeshPhongMaterial ? (t(e, n), function(e, t) {
				e.specular.value.copy(t.specular), e.shininess.value = Math.max(t.shininess, 1e-4), t.emissiveMap && (e.emissiveMap.value = t.emissiveMap);
				t.bumpMap && (e.bumpMap.value = t.bumpMap, e.bumpScale.value = t.bumpScale, 1 === t.side && (e.bumpScale.value *= -1));
				t.normalMap && (e.normalMap.value = t.normalMap, e.normalScale.value.copy(t.normalScale), 1 === t.side && e.normalScale.value.negate());
				t.displacementMap && (e.displacementMap.value = t.displacementMap, e.displacementScale.value = t.displacementScale, e.displacementBias.value = t.displacementBias)
			}(e, n)) : n.isMeshStandardMaterial ? (t(e, n), n.isMeshPhysicalMaterial ? function(e, t) {
				i(e, t), e.reflectivity.value = t.reflectivity, e.clearcoat.value = t.clearcoat, e.clearcoatRoughness.value = t.clearcoatRoughness, t.sheen && e.sheen.value.copy(t.sheen);
				t.clearcoatMap && (e.clearcoatMap.value = t.clearcoatMap);
				t.clearcoatRoughnessMap && (e.clearcoatRoughnessMap.value = t.clearcoatRoughnessMap);
				t.clearcoatNormalMap && (e.clearcoatNormalScale.value.copy(t.clearcoatNormalScale), e.clearcoatNormalMap.value = t.clearcoatNormalMap, 1 === t.side && e.clearcoatNormalScale.value.negate());
				e.transmission.value = t.transmission, t.transmissionMap && (e.transmissionMap.value = t.transmissionMap)
			}(e, n) : i(e, n)) : n.isMeshMatcapMaterial ? (t(e, n), function(e, t) {
				t.matcap && (e.matcap.value = t.matcap);
				t.bumpMap && (e.bumpMap.value = t.bumpMap, e.bumpScale.value = t.bumpScale, 1 === t.side && (e.bumpScale.value *= -1));
				t.normalMap && (e.normalMap.value = t.normalMap, e.normalScale.value.copy(t.normalScale), 1 === t.side && e.normalScale.value.negate());
				t.displacementMap && (e.displacementMap.value = t.displacementMap, e.displacementScale.value = t.displacementScale, e.displacementBias.value = t.displacementBias)
			}(e, n)) : n.isMeshDepthMaterial ? (t(e, n), function(e, t) {
				t.displacementMap && (e.displacementMap.value = t.displacementMap, e.displacementScale.value = t.displacementScale, e.displacementBias.value = t.displacementBias)
			}(e, n)) : n.isMeshDistanceMaterial ? (t(e, n), function(e, t) {
				t.displacementMap && (e.displacementMap.value = t.displacementMap, e.displacementScale.value = t.displacementScale, e.displacementBias.value = t.displacementBias);
				e.referencePosition.value.copy(t.referencePosition), e.nearDistance.value = t.nearDistance, e.farDistance.value = t.farDistance
			}(e, n)) : n.isMeshNormalMaterial ? (t(e, n), function(e, t) {
				t.bumpMap && (e.bumpMap.value = t.bumpMap, e.bumpScale.value = t.bumpScale, 1 === t.side && (e.bumpScale.value *= -1));
				t.normalMap && (e.normalMap.value = t.normalMap, e.normalScale.value.copy(t.normalScale), 1 === t.side && e.normalScale.value.negate());
				t.displacementMap && (e.displacementMap.value = t.displacementMap, e.displacementScale.value = t.displacementScale, e.displacementBias.value = t.displacementBias)
			}(e, n)) : n.isLineBasicMaterial ? (function(e, t) {
				e.diffuse.value.copy(t.color), e.opacity.value = t.opacity
			}(e, n), n.isLineDashedMaterial && function(e, t) {
				e.dashSize.value = t.dashSize, e.totalSize.value = t.dashSize + t.gapSize, e.scale.value = t.scale
			}(e, n)) : n.isPointsMaterial ? function(e, t, i, n) {
				e.diffuse.value.copy(t.color), e.opacity.value = t.opacity, e.size.value = t.size * i, e.scale.value = .5 * n, t.map && (e.map.value = t.map);
				t.alphaMap && (e.alphaMap.value = t.alphaMap);
				let r;
				t.map ? r = t.map : t.alphaMap && (r = t.alphaMap);
				void 0 !== r && (!0 === r.matrixAutoUpdate && r.updateMatrix(), e.uvTransform.value.copy(r.matrix))
			}(e, n, r, a) : n.isSpriteMaterial ? function(e, t) {
				e.diffuse.value.copy(t.color), e.opacity.value = t.opacity, e.rotation.value = t.rotation, t.map && (e.map.value = t.map);
				t.alphaMap && (e.alphaMap.value = t.alphaMap);
				let i;
				t.map ? i = t.map : t.alphaMap && (i = t.alphaMap);
				void 0 !== i && (!0 === i.matrixAutoUpdate && i.updateMatrix(), e.uvTransform.value.copy(i.matrix))
			}(e, n) : n.isShadowMaterial ? (e.color.value.copy(n.color), e.opacity.value = n.opacity) : n.isShaderMaterial && (n.uniformsNeedUpdate = !1)
		}
	}
}

function createCanvasElement() {
	const e = document.createElementNS("http://www.w3.org/1999/xhtml", "canvas");
	return e.style.display = "block", e
}

function WebGLRenderer(e) {
	const t = void 0 !== (e = e || {}).canvas ? e.canvas : createCanvasElement(),
		i = void 0 !== e.context ? e.context : null,
		n = void 0 !== e.alpha && e.alpha,
		r = void 0 === e.depth || e.depth,
		a = void 0 === e.stencil || e.stencil,
		s = void 0 !== e.antialias && e.antialias,
		o = void 0 === e.premultipliedAlpha || e.premultipliedAlpha,
		l = void 0 !== e.preserveDrawingBuffer && e.preserveDrawingBuffer,
		c = void 0 !== e.powerPreference ? e.powerPreference : "default",
		h = void 0 !== e.failIfMajorPerformanceCaveat && e.failIfMajorPerformanceCaveat;
	let u = null,
		d = null;
	const p = [],
		m = [];
	this.domElement = t, this.debug = {
		checkShaderErrors: !0
	}, this.autoClear = !0, this.autoClearColor = !0, this.autoClearDepth = !0, this.autoClearStencil = !0, this.sortObjects = !0, this.clippingPlanes = [], this.localClippingEnabled = !1, this.gammaFactor = 2, this.outputEncoding = 3e3, this.physicallyCorrectLights = !1, this.toneMapping = 0, this.toneMappingExposure = 1, this.maxMorphTargets = 8, this.maxMorphNormals = 4;
	const A = this;
	let g = !1,
		f = null,
		v = 0,
		y = 0,
		E = null,
		_ = null,
		b = -1,
		x = null;
	const w = new Vector4,
		C = new Vector4;
	let S = null,
		I = t.width,
		M = t.height,
		T = 1,
		B = null,
		L = null;
	const R = new Vector4(0, 0, I, M),
		D = new Vector4(0, 0, I, M);
	let P = !1;
	const Q = new Frustum;
	let F = !1,
		O = !1;
	const N = new Matrix4,
		k = new Vector3,
		U = {
			background: null,
			fog: null,
			environment: null,
			overrideMaterial: null,
			isScene: !0
		};

	function G() {
		return null === E ? T : 1
	}
	let H, V, $, z, W, q, j, Y, X, J, K, Z, ee, te, ie, ne, re, ae, se, oe, le, ce = i;

	function he(e, i) {
		for (let n = 0; n < e.length; n++) {
			const r = e[n],
				a = t.getContext(r, i);
			if (null !== a) return a
		}
		return null
	}
	try {
		const e = {
			alpha: n,
			depth: r,
			stencil: a,
			antialias: s,
			premultipliedAlpha: o,
			preserveDrawingBuffer: l,
			powerPreference: c,
			failIfMajorPerformanceCaveat: h
		};
		if (t.addEventListener("webglcontextlost", me, !1), t.addEventListener("webglcontextrestored", Ae, !1), null === ce) {
			const t = ["webgl2", "webgl", "experimental-webgl"];
			if (!0 === A.isWebGL1Renderer && t.shift(), ce = he(t, e), null === ce) throw he(t) ? new Error("Error creating WebGL context with your selected attributes.") : new Error("Error creating WebGL context.")
		}
		void 0 === ce.getShaderPrecisionFormat && (ce.getShaderPrecisionFormat = function() {
			return {
				rangeMin: 1,
				rangeMax: 1,
				precision: 1
			}
		})
	} catch (e) {
		throw console.error("THREE.WebGLRenderer: " + e.message), e
	}

	function ue() {
		H = new WebGLExtensions(ce), V = new WebGLCapabilities(ce, H, e), H.init(V), oe = new WebGLUtils(ce, H, V), $ = new WebGLState(ce, H, V), $.scissor(C.copy(D).multiplyScalar(T).floor()), $.viewport(w.copy(R).multiplyScalar(T).floor()), z = new WebGLInfo(ce), W = new WebGLProperties, q = new WebGLTextures(ce, H, $, W, V, oe, z), j = new WebGLCubeMaps(A), Y = new WebGLAttributes(ce, V), le = new WebGLBindingStates(ce, H, Y, V), X = new WebGLGeometries(ce, Y, z, le), J = new WebGLObjects(ce, X, Y, z), re = new WebGLMorphtargets(ce), ie = new WebGLClipping(W), K = new WebGLPrograms(A, j, H, V, le, ie), Z = new WebGLMaterials(W), ee = new WebGLRenderLists(W), te = new WebGLRenderStates(H, V), ne = new WebGLBackground(A, j, $, J, o), ae = new WebGLBufferRenderer(ce, H, z, V), se = new WebGLIndexedBufferRenderer(ce, H, z, V), z.programs = K.programs, A.capabilities = V, A.extensions = H, A.properties = W, A.renderLists = ee, A.state = $, A.info = z
	}
	ue();
	const de = new WebXRManager(A, ce);
	this.xr = de;
	const pe = new WebGLShadowMap(A, J, V.maxTextureSize);

	function me(e) {
		e.preventDefault(), console.log("THREE.WebGLRenderer: Context Lost."), g = !0
	}

	function Ae() {
		console.log("THREE.WebGLRenderer: Context Restored."), g = !1, ue()
	}

	function ge(e) {
		const t = e.target;
		t.removeEventListener("dispose", ge),
			function(e) {
				fe(e), W.remove(e)
			}(t)
	}

	function fe(e) {
		const t = W.get(e).program;
		void 0 !== t && K.releaseProgram(t)
	}
	this.shadowMap = pe, this.getContext = function() {
		return ce
	}, this.getContextAttributes = function() {
		return ce.getContextAttributes()
	}, this.forceContextLoss = function() {
		const e = H.get("WEBGL_lose_context");
		e && e.loseContext()
	}, this.forceContextRestore = function() {
		const e = H.get("WEBGL_lose_context");
		e && e.restoreContext()
	}, this.getPixelRatio = function() {
		return T
	}, this.setPixelRatio = function(e) {
		void 0 !== e && (T = e, this.setSize(I, M, !1))
	}, this.getSize = function(e) {
		return void 0 === e && (console.warn("WebGLRenderer: .getsize() now requires a Vector2 as an argument"), e = new Vector2), e.set(I, M)
	}, this.setSize = function(e, i, n) {
		de.isPresenting ? console.warn("THREE.WebGLRenderer: Can't change size while VR device is presenting.") : (I = e, M = i, t.width = Math.floor(e * T), t.height = Math.floor(i * T), !1 !== n && (t.style.width = e + "px", t.style.height = i + "px"), this.setViewport(0, 0, e, i))
	}, this.getDrawingBufferSize = function(e) {
		return void 0 === e && (console.warn("WebGLRenderer: .getdrawingBufferSize() now requires a Vector2 as an argument"), e = new Vector2), e.set(I * T, M * T).floor()
	}, this.setDrawingBufferSize = function(e, i, n) {
		I = e, M = i, T = n, t.width = Math.floor(e * n), t.height = Math.floor(i * n), this.setViewport(0, 0, e, i)
	}, this.getCurrentViewport = function(e) {
		return void 0 === e && (console.warn("WebGLRenderer: .getCurrentViewport() now requires a Vector4 as an argument"), e = new Vector4), e.copy(w)
	}, this.getViewport = function(e) {
		return e.copy(R)
	}, this.setViewport = function(e, t, i, n) {
		e.isVector4 ? R.set(e.x, e.y, e.z, e.w) : R.set(e, t, i, n), $.viewport(w.copy(R).multiplyScalar(T).floor())
	}, this.getScissor = function(e) {
		return e.copy(D)
	}, this.setScissor = function(e, t, i, n) {
		e.isVector4 ? D.set(e.x, e.y, e.z, e.w) : D.set(e, t, i, n), $.scissor(C.copy(D).multiplyScalar(T).floor())
	}, this.getScissorTest = function() {
		return P
	}, this.setScissorTest = function(e) {
		$.setScissorTest(P = e)
	}, this.setOpaqueSort = function(e) {
		B = e
	}, this.setTransparentSort = function(e) {
		L = e
	}, this.getClearColor = function(e) {
		return void 0 === e && (console.warn("WebGLRenderer: .getClearColor() now requires a Color as an argument"), e = new Color), e.copy(ne.getClearColor())
	}, this.setClearColor = function() {
		ne.setClearColor.apply(ne, arguments)
	}, this.getClearAlpha = function() {
		return ne.getClearAlpha()
	}, this.setClearAlpha = function() {
		ne.setClearAlpha.apply(ne, arguments)
	}, this.clear = function(e, t, i) {
		let n = 0;
		(void 0 === e || e) && (n |= 16384), (void 0 === t || t) && (n |= 256), (void 0 === i || i) && (n |= 1024), ce.clear(n)
	}, this.clearColor = function() {
		this.clear(!0, !1, !1)
	}, this.clearDepth = function() {
		this.clear(!1, !0, !1)
	}, this.clearStencil = function() {
		this.clear(!1, !1, !0)
	}, this.dispose = function() {
		t.removeEventListener("webglcontextlost", me, !1), t.removeEventListener("webglcontextrestored", Ae, !1), ee.dispose(), te.dispose(), W.dispose(), j.dispose(), J.dispose(), le.dispose(), de.dispose(), ye.stop()
	}, this.renderBufferImmediate = function(e, t) {
		le.initAttributes();
		const i = W.get(e);
		e.hasPositions && !i.position && (i.position = ce.createBuffer()), e.hasNormals && !i.normal && (i.normal = ce.createBuffer()), e.hasUvs && !i.uv && (i.uv = ce.createBuffer()), e.hasColors && !i.color && (i.color = ce.createBuffer());
		const n = t.getAttributes();
		e.hasPositions && (ce.bindBuffer(34962, i.position), ce.bufferData(34962, e.positionArray, 35048), le.enableAttribute(n.position), ce.vertexAttribPointer(n.position, 3, 5126, !1, 0, 0)), e.hasNormals && (ce.bindBuffer(34962, i.normal), ce.bufferData(34962, e.normalArray, 35048), le.enableAttribute(n.normal), ce.vertexAttribPointer(n.normal, 3, 5126, !1, 0, 0)), e.hasUvs && (ce.bindBuffer(34962, i.uv), ce.bufferData(34962, e.uvArray, 35048), le.enableAttribute(n.uv), ce.vertexAttribPointer(n.uv, 2, 5126, !1, 0, 0)), e.hasColors && (ce.bindBuffer(34962, i.color), ce.bufferData(34962, e.colorArray, 35048), le.enableAttribute(n.color), ce.vertexAttribPointer(n.color, 3, 5126, !1, 0, 0)), le.disableUnusedAttributes(), ce.drawArrays(4, 0, e.count), e.count = 0
	}, this.renderBufferDirect = function(e, t, i, n, r, a) {
		null === t && (t = U);
		const s = r.isMesh && r.matrixWorld.determinant() < 0,
			o = we(e, t, n, r);
		$.setMaterial(n, s);
		let l = i.index;
		const c = i.attributes.position;
		if (null === l) {
			if (void 0 === c || 0 === c.count) return
		} else if (0 === l.count) return;
		let h, u = 1;
		!0 === n.wireframe && (l = X.getWireframeAttribute(i), u = 2), (n.morphTargets || n.morphNormals) && re.update(r, i, n, o), le.setup(r, n, o, i, l);
		let d = ae;
		null !== l && (h = Y.get(l), d = se, d.setIndex(h));
		const p = null !== l ? l.count : c.count,
			m = i.drawRange.start * u,
			A = i.drawRange.count * u,
			g = null !== a ? a.start * u : 0,
			f = null !== a ? a.count * u : 1 / 0,
			v = Math.max(m, g),
			y = Math.min(p, m + A, g + f) - 1,
			E = Math.max(0, y - v + 1);
		if (0 !== E) {
			if (r.isMesh) !0 === n.wireframe ? ($.setLineWidth(n.wireframeLinewidth * G()), d.setMode(1)) : d.setMode(4);
			else if (r.isLine) {
				let e = n.linewidth;
				void 0 === e && (e = 1), $.setLineWidth(e * G()), r.isLineSegments ? d.setMode(1) : r.isLineLoop ? d.setMode(2) : d.setMode(3)
			} else r.isPoints ? d.setMode(0) : r.isSprite && d.setMode(4);
			if (r.isInstancedMesh) d.renderInstances(v, E, r.count);
			else if (i.isInstancedBufferGeometry) {
				const e = Math.min(i.instanceCount, i._maxInstanceCount);
				d.renderInstances(v, E, e)
			} else d.render(v, E)
		}
	}, this.compile = function(e, t) {
		d = te.get(e), d.init(), e.traverseVisible((function(e) {
			e.isLight && e.layers.test(t.layers) && (d.pushLight(e), e.castShadow && d.pushShadow(e))
		})), d.setupLights();
		const i = new WeakMap;
		e.traverse((function(t) {
			const n = t.material;
			if (n)
				if (Array.isArray(n))
					for (let r = 0; r < n.length; r++) {
						const a = n[r];
						!1 === i.has(a) && (xe(a, e, t), i.set(a))
					} else !1 === i.has(n) && (xe(n, e, t), i.set(n))
		}))
	};
	let ve = null;
	const ye = new WebGLAnimation;

	function Ee(e, t, i, n) {
		if (!1 === e.visible) return;
		if (e.layers.test(t.layers))
			if (e.isGroup) i = e.renderOrder;
			else if (e.isLOD) !0 === e.autoUpdate && e.update(t);
		else if (e.isLight) d.pushLight(e), e.castShadow && d.pushShadow(e);
		else if (e.isSprite) {
			if (!e.frustumCulled || Q.intersectsSprite(e)) {
				n && k.setFromMatrixPosition(e.matrixWorld).applyMatrix4(N);
				const t = J.update(e),
					r = e.material;
				r.visible && u.push(e, t, r, i, k.z, null)
			}
		} else if (e.isImmediateRenderObject) n && k.setFromMatrixPosition(e.matrixWorld).applyMatrix4(N), u.push(e, null, e.material, i, k.z, null);
		else if ((e.isMesh || e.isLine || e.isPoints) && (e.isSkinnedMesh && e.skeleton.frame !== z.render.frame && (e.skeleton.update(), e.skeleton.frame = z.render.frame), !e.frustumCulled || Q.intersectsObject(e))) {
			n && k.setFromMatrixPosition(e.matrixWorld).applyMatrix4(N);
			const t = J.update(e),
				r = e.material;
			if (Array.isArray(r)) {
				const n = t.groups;
				for (let a = 0, s = n.length; a < s; a++) {
					const s = n[a],
						o = r[s.materialIndex];
					o && o.visible && u.push(e, t, o, i, k.z, s)
				}
			} else r.visible && u.push(e, t, r, i, k.z, null)
		}
		const r = e.children;
		for (let e = 0, a = r.length; e < a; e++) Ee(r[e], t, i, n)
	}

	function _e(e, t, i) {
		const n = !0 === t.isScene ? t.overrideMaterial : null;
		for (let r = 0, a = e.length; r < a; r++) {
			const a = e[r],
				s = a.object,
				o = a.geometry,
				l = null === n ? a.material : n,
				c = a.group;
			if (i.isArrayCamera) {
				const e = i.cameras;
				for (let i = 0, n = e.length; i < n; i++) {
					const n = e[i];
					s.layers.test(n.layers) && ($.viewport(w.copy(n.viewport)), d.setupLightsView(n), be(s, t, n, o, l, c))
				}
			} else be(s, t, i, o, l, c)
		}
	}

	function be(e, t, i, n, r, a) {
		if (e.onBeforeRender(A, t, i, n, r, a), e.modelViewMatrix.multiplyMatrices(i.matrixWorldInverse, e.matrixWorld), e.normalMatrix.getNormalMatrix(e.modelViewMatrix), e.isImmediateRenderObject) {
			const n = we(i, t, r, e);
			$.setMaterial(r), le.reset(),
				function(e, t) {
					e.render((function(e) {
						A.renderBufferImmediate(e, t)
					}))
				}(e, n)
		} else A.renderBufferDirect(i, t, n, r, e, a);
		e.onAfterRender(A, t, i, n, r, a)
	}

	function xe(e, t, i) {
		!0 !== t.isScene && (t = U);
		const n = W.get(e),
			r = d.state.lights,
			a = d.state.shadowsArray,
			s = r.state.version,
			o = K.getParameters(e, r.state, a, t, i),
			l = K.getProgramCacheKey(o);
		let c = n.program,
			h = !0;
		if (n.environment = e.isMeshStandardMaterial ? t.environment : null, n.fog = t.fog, n.envMap = j.get(e.envMap || n.environment), void 0 === c) e.addEventListener("dispose", ge);
		else if (c.cacheKey !== l) fe(e);
		else if (n.lightsStateVersion !== s) h = !1;
		else {
			if (void 0 !== o.shaderID) return;
			h = !1
		}
		h && (o.uniforms = K.getUniforms(e), e.onBeforeCompile(o, A), c = K.acquireProgram(o, l), n.program = c, n.uniforms = o.uniforms, n.outputEncoding = o.outputEncoding);
		const u = n.uniforms;
		(e.isShaderMaterial || e.isRawShaderMaterial) && !0 !== e.clipping || (n.numClippingPlanes = ie.numPlanes, n.numIntersection = ie.numIntersection, u.clippingPlanes = ie.uniform), n.needsLights = function(e) {
			return e.isMeshLambertMaterial || e.isMeshToonMaterial || e.isMeshPhongMaterial || e.isMeshStandardMaterial || e.isShadowMaterial || e.isShaderMaterial && !0 === e.lights
		}(e), n.lightsStateVersion = s, n.needsLights && (u.ambientLightColor.value = r.state.ambient, u.lightProbe.value = r.state.probe, u.directionalLights.value = r.state.directional, u.directionalLightShadows.value = r.state.directionalShadow, u.spotLights.value = r.state.spot, u.spotLightShadows.value = r.state.spotShadow, u.rectAreaLights.value = r.state.rectArea, u.ltc_1.value = r.state.rectAreaLTC1, u.ltc_2.value = r.state.rectAreaLTC2, u.pointLights.value = r.state.point, u.pointLightShadows.value = r.state.pointShadow, u.hemisphereLights.value = r.state.hemi, u.directionalShadowMap.value = r.state.directionalShadowMap, u.directionalShadowMatrix.value = r.state.directionalShadowMatrix, u.spotShadowMap.value = r.state.spotShadowMap, u.spotShadowMatrix.value = r.state.spotShadowMatrix, u.pointShadowMap.value = r.state.pointShadowMap, u.pointShadowMatrix.value = r.state.pointShadowMatrix);
		const p = n.program.getUniforms(),
			m = WebGLUniforms.seqWithValue(p.seq, u);
		n.uniformsList = m
	}

	function we(e, t, i, n) {
		!0 !== t.isScene && (t = U), q.resetTextureUnits();
		const r = t.fog,
			a = i.isMeshStandardMaterial ? t.environment : null,
			s = null === E ? A.outputEncoding : E.texture.encoding,
			o = j.get(i.envMap || a),
			l = W.get(i),
			c = d.state.lights;
		if (!0 === F && (!0 === O || e !== x)) {
			const t = e === x && i.id === b;
			ie.setState(i, e, t)
		}
		i.version === l.__version ? i.fog && l.fog !== r || l.environment !== a || l.needsLights && l.lightsStateVersion !== c.state.version ? xe(i, t, n) : void 0 === l.numClippingPlanes || l.numClippingPlanes === ie.numPlanes && l.numIntersection === ie.numIntersection ? (l.outputEncoding !== s || l.envMap !== o) && xe(i, t, n) : xe(i, t, n) : (xe(i, t, n), l.__version = i.version);
		let h = !1,
			u = !1,
			p = !1;
		const m = l.program,
			g = m.getUniforms(),
			f = l.uniforms;
		if ($.useProgram(m.program) && (h = !0, u = !0, p = !0), i.id !== b && (b = i.id, u = !0), h || x !== e) {
			if (g.setValue(ce, "projectionMatrix", e.projectionMatrix), V.logarithmicDepthBuffer && g.setValue(ce, "logDepthBufFC", 2 / (Math.log(e.far + 1) / Math.LN2)), x !== e && (x = e, u = !0, p = !0), i.isShaderMaterial || i.isMeshPhongMaterial || i.isMeshToonMaterial || i.isMeshStandardMaterial || i.envMap) {
				const t = g.map.cameraPosition;
				void 0 !== t && t.setValue(ce, k.setFromMatrixPosition(e.matrixWorld))
			}(i.isMeshPhongMaterial || i.isMeshToonMaterial || i.isMeshLambertMaterial || i.isMeshBasicMaterial || i.isMeshStandardMaterial || i.isShaderMaterial) && g.setValue(ce, "isOrthographic", !0 === e.isOrthographicCamera), (i.isMeshPhongMaterial || i.isMeshToonMaterial || i.isMeshLambertMaterial || i.isMeshBasicMaterial || i.isMeshStandardMaterial || i.isShaderMaterial || i.isShadowMaterial || i.skinning) && g.setValue(ce, "viewMatrix", e.matrixWorldInverse)
		}
		if (i.skinning) {
			g.setOptional(ce, n, "bindMatrix"), g.setOptional(ce, n, "bindMatrixInverse");
			const e = n.skeleton;
			if (e) {
				const t = e.bones;
				if (V.floatVertexTextures) {
					if (null === e.boneTexture) {
						let i = Math.sqrt(4 * t.length);
						i = MathUtils.ceilPowerOfTwo(i), i = Math.max(i, 4);
						const n = new Float32Array(i * i * 4);
						n.set(e.boneMatrices);
						const r = new DataTexture(n, i, i, 1023, 1015);
						e.boneMatrices = n, e.boneTexture = r, e.boneTextureSize = i
					}
					g.setValue(ce, "boneTexture", e.boneTexture, q), g.setValue(ce, "boneTextureSize", e.boneTextureSize)
				} else g.setOptional(ce, e, "boneMatrices")
			}
		}
		var v, y;
		return (u || l.receiveShadow !== n.receiveShadow) && (l.receiveShadow = n.receiveShadow, g.setValue(ce, "receiveShadow", n.receiveShadow)), u && (g.setValue(ce, "toneMappingExposure", A.toneMappingExposure), l.needsLights && (y = p, (v = f).ambientLightColor.needsUpdate = y, v.lightProbe.needsUpdate = y, v.directionalLights.needsUpdate = y, v.directionalLightShadows.needsUpdate = y, v.pointLights.needsUpdate = y, v.pointLightShadows.needsUpdate = y, v.spotLights.needsUpdate = y, v.spotLightShadows.needsUpdate = y, v.rectAreaLights.needsUpdate = y, v.hemisphereLights.needsUpdate = y), r && i.fog && Z.refreshFogUniforms(f, r), Z.refreshMaterialUniforms(f, i, T, M), WebGLUniforms.upload(ce, l.uniformsList, f, q)), i.isShaderMaterial && !0 === i.uniformsNeedUpdate && (WebGLUniforms.upload(ce, l.uniformsList, f, q), i.uniformsNeedUpdate = !1), i.isSpriteMaterial && g.setValue(ce, "center", n.center), g.setValue(ce, "modelViewMatrix", n.modelViewMatrix), g.setValue(ce, "normalMatrix", n.normalMatrix), g.setValue(ce, "modelMatrix", n.matrixWorld), m
	}
	ye.setAnimationLoop((function(e) {
		de.isPresenting || ve && ve(e)
	})), "undefined" != typeof window && ye.setContext(window), this.setAnimationLoop = function(e) {
		ve = e, de.setAnimationLoop(e), null === e ? ye.stop() : ye.start()
	}, this.render = function(e, t) {
		let i, n;
		if (void 0 !== arguments[2] && (console.warn("THREE.WebGLRenderer.render(): the renderTarget argument has been removed. Use .setRenderTarget() instead."), i = arguments[2]), void 0 !== arguments[3] && (console.warn("THREE.WebGLRenderer.render(): the forceClear argument has been removed. Use .clear() instead."), n = arguments[3]), void 0 !== t && !0 !== t.isCamera) return void console.error("THREE.WebGLRenderer.render: camera is not an instance of THREE.Camera.");
		if (!0 === g) return;
		le.resetDefaultState(), b = -1, x = null, !0 === e.autoUpdate && e.updateMatrixWorld(), null === t.parent && t.updateMatrixWorld(), !0 === de.enabled && !0 === de.isPresenting && (t = de.getCamera(t)), !0 === e.isScene && e.onBeforeRender(A, e, t, i || E), d = te.get(e, m.length), d.init(), m.push(d), N.multiplyMatrices(t.projectionMatrix, t.matrixWorldInverse), Q.setFromProjectionMatrix(N), O = this.localClippingEnabled, F = ie.init(this.clippingPlanes, O, t), u = ee.get(e, p.length), u.init(), p.push(u), Ee(e, t, 0, A.sortObjects), u.finish(), !0 === A.sortObjects && u.sort(B, L), !0 === F && ie.beginShadows();
		const r = d.state.shadowsArray;
		pe.render(r, e, t), d.setupLights(), d.setupLightsView(t), !0 === F && ie.endShadows(), !0 === this.info.autoReset && this.info.reset(), void 0 !== i && this.setRenderTarget(i), ne.render(u, e, t, n);
		const a = u.opaque,
			s = u.transparent;
		a.length > 0 && _e(a, e, t), s.length > 0 && _e(s, e, t), !0 === e.isScene && e.onAfterRender(A, e, t), null !== E && (q.updateRenderTargetMipmap(E), q.updateMultisampleRenderTarget(E)), $.buffers.depth.setTest(!0), $.buffers.depth.setMask(!0), $.buffers.color.setMask(!0), $.setPolygonOffset(!1), m.pop(), d = m.length > 0 ? m[m.length - 1] : null, p.pop(), u = p.length > 0 ? p[p.length - 1] : null
	}, this.setFramebuffer = function(e) {
		f !== e && null === E && ce.bindFramebuffer(36160, e), f = e
	}, this.getActiveCubeFace = function() {
		return v
	}, this.getActiveMipmapLevel = function() {
		return y
	}, this.getRenderTarget = function() {
		return E
	}, this.setRenderTarget = function(e, t = 0, i = 0) {
		E = e, v = t, y = i, e && void 0 === W.get(e).__webglFramebuffer && q.setupRenderTarget(e);
		let n = f,
			r = !1,
			a = !1;
		if (e) {
			const i = e.texture;
			(i.isDataTexture3D || i.isDataTexture2DArray) && (a = !0);
			const s = W.get(e).__webglFramebuffer;
			e.isWebGLCubeRenderTarget ? (n = s[t], r = !0) : n = e.isWebGLMultisampleRenderTarget ? W.get(e).__webglMultisampledFramebuffer : s, w.copy(e.viewport), C.copy(e.scissor), S = e.scissorTest
		} else w.copy(R).multiplyScalar(T).floor(), C.copy(D).multiplyScalar(T).floor(), S = P;
		if (_ !== n && (ce.bindFramebuffer(36160, n), _ = n), $.viewport(w), $.scissor(C), $.setScissorTest(S), r) {
			const n = W.get(e.texture);
			ce.framebufferTexture2D(36160, 36064, 34069 + t, n.__webglTexture, i)
		} else if (a) {
			const n = W.get(e.texture),
				r = t || 0;
			ce.framebufferTextureLayer(36160, 36064, n.__webglTexture, i || 0, r)
		}
	}, this.readRenderTargetPixels = function(e, t, i, n, r, a, s) {
		if (!e || !e.isWebGLRenderTarget) return void console.error("THREE.WebGLRenderer.readRenderTargetPixels: renderTarget is not THREE.WebGLRenderTarget.");
		let o = W.get(e).__webglFramebuffer;
		if (e.isWebGLCubeRenderTarget && void 0 !== s && (o = o[s]), o) {
			let s = !1;
			o !== _ && (ce.bindFramebuffer(36160, o), s = !0);
			try {
				const o = e.texture,
					l = o.format,
					c = o.type;
				if (1023 !== l && oe.convert(l) !== ce.getParameter(35739)) return void console.error("THREE.WebGLRenderer.readRenderTargetPixels: renderTarget is not in RGBA or implementation defined format.");
				const h = 1016 === c && (H.has("EXT_color_buffer_half_float") || V.isWebGL2 && H.has("EXT_color_buffer_float"));
				if (!(1009 === c || oe.convert(c) === ce.getParameter(35738) || 1015 === c && (V.isWebGL2 || H.has("OES_texture_float") || H.has("WEBGL_color_buffer_float")) || h)) return void console.error("THREE.WebGLRenderer.readRenderTargetPixels: renderTarget is not in UnsignedByteType or implementation defined type.");
				36053 === ce.checkFramebufferStatus(36160) ? t >= 0 && t <= e.width - n && i >= 0 && i <= e.height - r && ce.readPixels(t, i, n, r, oe.convert(l), oe.convert(c), a) : console.error("THREE.WebGLRenderer.readRenderTargetPixels: readPixels from renderTarget failed. Framebuffer not complete.")
			} finally {
				s && ce.bindFramebuffer(36160, _)
			}
		}
	}, this.copyFramebufferToTexture = function(e, t, i = 0) {
		const n = Math.pow(2, -i),
			r = Math.floor(t.image.width * n),
			a = Math.floor(t.image.height * n),
			s = oe.convert(t.format);
		q.setTexture2D(t, 0), ce.copyTexImage2D(3553, i, s, e.x, e.y, r, a, 0), $.unbindTexture()
	}, this.copyTextureToTexture = function(e, t, i, n = 0) {
		const r = t.image.width,
			a = t.image.height,
			s = oe.convert(i.format),
			o = oe.convert(i.type);
		q.setTexture2D(i, 0), ce.pixelStorei(37440, i.flipY), ce.pixelStorei(37441, i.premultiplyAlpha), ce.pixelStorei(3317, i.unpackAlignment), t.isDataTexture ? ce.texSubImage2D(3553, n, e.x, e.y, r, a, s, o, t.image.data) : t.isCompressedTexture ? ce.compressedTexSubImage2D(3553, n, e.x, e.y, t.mipmaps[0].width, t.mipmaps[0].height, s, t.mipmaps[0].data) : ce.texSubImage2D(3553, n, e.x, e.y, s, o, t.image), 0 === n && i.generateMipmaps && ce.generateMipmap(3553), $.unbindTexture()
	}, this.copyTextureToTexture3D = function(e, t, i, n, r = 0) {
		if (A.isWebGL1Renderer) return void console.warn("THREE.WebGLRenderer.copyTextureToTexture3D: can only be used with WebGL2.");
		const {
			width: a,
			height: s,
			data: o
		} = i.image, l = oe.convert(n.format), c = oe.convert(n.type);
		let h;
		if (n.isDataTexture3D) q.setTexture3D(n, 0), h = 32879;
		else {
			if (!n.isDataTexture2DArray) return void console.warn("THREE.WebGLRenderer.copyTextureToTexture3D: only supports THREE.DataTexture3D and THREE.DataTexture2DArray.");
			q.setTexture2DArray(n, 0), h = 35866
		}
		ce.pixelStorei(37440, n.flipY), ce.pixelStorei(37441, n.premultiplyAlpha), ce.pixelStorei(3317, n.unpackAlignment);
		const u = ce.getParameter(3314),
			d = ce.getParameter(32878),
			p = ce.getParameter(3316),
			m = ce.getParameter(3315),
			g = ce.getParameter(32877);
		ce.pixelStorei(3314, a), ce.pixelStorei(32878, s), ce.pixelStorei(3316, e.min.x), ce.pixelStorei(3315, e.min.y), ce.pixelStorei(32877, e.min.z), ce.texSubImage3D(h, r, t.x, t.y, t.z, e.max.x - e.min.x + 1, e.max.y - e.min.y + 1, e.max.z - e.min.z + 1, l, c, o), ce.pixelStorei(3314, u), ce.pixelStorei(32878, d), ce.pixelStorei(3316, p), ce.pixelStorei(3315, m), ce.pixelStorei(32877, g), 0 === r && n.generateMipmaps && ce.generateMipmap(h), $.unbindTexture()
	}, this.initTexture = function(e) {
		q.setTexture2D(e, 0), $.unbindTexture()
	}, this.resetState = function() {
		$.reset(), le.reset()
	}, "undefined" != typeof __THREE_DEVTOOLS__ && __THREE_DEVTOOLS__.dispatchEvent(new CustomEvent("observe", {
		detail: this
	}))
}
Group.prototype.isGroup = !0, Object.assign(WebXRController.prototype, {
	constructor: WebXRController,
	getHandSpace: function() {
		return null === this._hand && (this._hand = new Group, this._hand.matrixAutoUpdate = !1, this._hand.visible = !1, this._hand.joints = {}, this._hand.inputState = {
			pinching: !1
		}), this._hand
	},
	getTargetRaySpace: function() {
		return null === this._targetRay && (this._targetRay = new Group, this._targetRay.matrixAutoUpdate = !1, this._targetRay.visible = !1), this._targetRay
	},
	getGripSpace: function() {
		return null === this._grip && (this._grip = new Group, this._grip.matrixAutoUpdate = !1, this._grip.visible = !1), this._grip
	},
	dispatchEvent: function(e) {
		return null !== this._targetRay && this._targetRay.dispatchEvent(e), null !== this._grip && this._grip.dispatchEvent(e), null !== this._hand && this._hand.dispatchEvent(e), this
	},
	disconnect: function(e) {
		return this.dispatchEvent({
			type: "disconnected",
			data: e
		}), null !== this._targetRay && (this._targetRay.visible = !1), null !== this._grip && (this._grip.visible = !1), null !== this._hand && (this._hand.visible = !1), this
	},
	update: function(e, t, i) {
		let n = null,
			r = null,
			a = null;
		const s = this._targetRay,
			o = this._grip,
			l = this._hand;
		if (e && "visible-blurred" !== t.session.visibilityState)
			if (l && e.hand) {
				a = !0;
				for (const n of e.hand.values()) {
					const e = t.getJointPose(n, i);
					if (void 0 === l.joints[n.jointName]) {
						const e = new Group;
						e.matrixAutoUpdate = !1, e.visible = !1, l.joints[n.jointName] = e, l.add(e)
					}
					const r = l.joints[n.jointName];
					null !== e && (r.matrix.fromArray(e.transform.matrix), r.matrix.decompose(r.position, r.rotation, r.scale), r.jointRadius = e.radius), r.visible = null !== e
				}
				const n = l.joints["index-finger-tip"],
					r = l.joints["thumb-tip"],
					s = n.position.distanceTo(r.position),
					o = .02,
					c = .005;
				l.inputState.pinching && s > o + c ? (l.inputState.pinching = !1, this.dispatchEvent({
					type: "pinchend",
					handedness: e.handedness,
					target: this
				})) : !l.inputState.pinching && s <= o - c && (l.inputState.pinching = !0, this.dispatchEvent({
					type: "pinchstart",
					handedness: e.handedness,
					target: this
				}))
			} else null !== s && (n = t.getPose(e.targetRaySpace, i), null !== n && (s.matrix.fromArray(n.transform.matrix), s.matrix.decompose(s.position, s.rotation, s.scale))), null !== o && e.gripSpace && (r = t.getPose(e.gripSpace, i), null !== r && (o.matrix.fromArray(r.transform.matrix), o.matrix.decompose(o.position, o.rotation, o.scale)));
		return null !== s && (s.visible = null !== n), null !== o && (o.visible = null !== r), null !== l && (l.visible = null !== a), this
	}
}), Object.assign(WebXRManager.prototype, EventDispatcher.prototype);
class WebGL1Renderer extends WebGLRenderer {}
WebGL1Renderer.prototype.isWebGL1Renderer = !0;
class Scene extends Object3D {
	constructor() {
		super(), this.type = "Scene", this.background = null, this.environment = null, this.fog = null, this.overrideMaterial = null, this.autoUpdate = !0, "undefined" != typeof __THREE_DEVTOOLS__ && __THREE_DEVTOOLS__.dispatchEvent(new CustomEvent("observe", {
			detail: this
		}))
	}
	copy(e, t) {
		return super.copy(e, t), null !== e.background && (this.background = e.background.clone()), null !== e.environment && (this.environment = e.environment.clone()), null !== e.fog && (this.fog = e.fog.clone()), null !== e.overrideMaterial && (this.overrideMaterial = e.overrideMaterial.clone()), this.autoUpdate = e.autoUpdate, this.matrixAutoUpdate = e.matrixAutoUpdate, this
	}
	toJSON(e) {
		const t = super.toJSON(e);
		return null !== this.background && (t.object.background = this.background.toJSON(e)), null !== this.environment && (t.object.environment = this.environment.toJSON(e)), null !== this.fog && (t.object.fog = this.fog.toJSON()), t
	}
}

function InterleavedBuffer(e, t) {
	this.array = e, this.stride = t, this.count = void 0 !== e ? e.length / t : 0, this.usage = 35044, this.updateRange = {
		offset: 0,
		count: -1
	}, this.version = 0, this.uuid = MathUtils.generateUUID()
}
Scene.prototype.isScene = !0, Object.defineProperty(InterleavedBuffer.prototype, "needsUpdate", {
	set: function(e) {
		!0 === e && this.version++
	}
}), Object.assign(InterleavedBuffer.prototype, {
	isInterleavedBuffer: !0,
	onUploadCallback: function() {},
	setUsage: function(e) {
		return this.usage = e, this
	},
	copy: function(e) {
		return this.array = new e.array.constructor(e.array), this.count = e.count, this.stride = e.stride, this.usage = e.usage, this
	},
	copyAt: function(e, t, i) {
		e *= this.stride, i *= t.stride;
		for (let n = 0, r = this.stride; n < r; n++) this.array[e + n] = t.array[i + n];
		return this
	},
	set: function(e, t = 0) {
		return this.array.set(e, t), this
	},
	clone: function(e) {
		void 0 === e.arrayBuffers && (e.arrayBuffers = {}), void 0 === this.array.buffer._uuid && (this.array.buffer._uuid = MathUtils.generateUUID()), void 0 === e.arrayBuffers[this.array.buffer._uuid] && (e.arrayBuffers[this.array.buffer._uuid] = this.array.slice(0).buffer);
		const t = new InterleavedBuffer(new this.array.constructor(e.arrayBuffers[this.array.buffer._uuid]), this.stride);
		return t.setUsage(this.usage), t
	},
	onUpload: function(e) {
		return this.onUploadCallback = e, this
	},
	toJSON: function(e) {
		return void 0 === e.arrayBuffers && (e.arrayBuffers = {}), void 0 === this.array.buffer._uuid && (this.array.buffer._uuid = MathUtils.generateUUID()), void 0 === e.arrayBuffers[this.array.buffer._uuid] && (e.arrayBuffers[this.array.buffer._uuid] = Array.prototype.slice.call(new Uint32Array(this.array.buffer))), {
			uuid: this.uuid,
			buffer: this.array.buffer._uuid,
			type: this.array.constructor.name,
			stride: this.stride
		}
	}
});
const _vector$6 = new Vector3;

function InterleavedBufferAttribute(e, t, i, n) {
	this.name = "", this.data = e, this.itemSize = t, this.offset = i, this.normalized = !0 === n
}
Object.defineProperties(InterleavedBufferAttribute.prototype, {
	count: {
		get: function() {
			return this.data.count
		}
	},
	array: {
		get: function() {
			return this.data.array
		}
	},
	needsUpdate: {
		set: function(e) {
			this.data.needsUpdate = e
		}
	}
}), Object.assign(InterleavedBufferAttribute.prototype, {
	isInterleavedBufferAttribute: !0,
	applyMatrix4: function(e) {
		for (let t = 0, i = this.data.count; t < i; t++) _vector$6.x = this.getX(t), _vector$6.y = this.getY(t), _vector$6.z = this.getZ(t), _vector$6.applyMatrix4(e), this.setXYZ(t, _vector$6.x, _vector$6.y, _vector$6.z);
		return this
	},
	setX: function(e, t) {
		return this.data.array[e * this.data.stride + this.offset] = t, this
	},
	setY: function(e, t) {
		return this.data.array[e * this.data.stride + this.offset + 1] = t, this
	},
	setZ: function(e, t) {
		return this.data.array[e * this.data.stride + this.offset + 2] = t, this
	},
	setW: function(e, t) {
		return this.data.array[e * this.data.stride + this.offset + 3] = t, this
	},
	getX: function(e) {
		return this.data.array[e * this.data.stride + this.offset]
	},
	getY: function(e) {
		return this.data.array[e * this.data.stride + this.offset + 1]
	},
	getZ: function(e) {
		return this.data.array[e * this.data.stride + this.offset + 2]
	},
	getW: function(e) {
		return this.data.array[e * this.data.stride + this.offset + 3]
	},
	setXY: function(e, t, i) {
		return e = e * this.data.stride + this.offset, this.data.array[e + 0] = t, this.data.array[e + 1] = i, this
	},
	setXYZ: function(e, t, i, n) {
		return e = e * this.data.stride + this.offset, this.data.array[e + 0] = t, this.data.array[e + 1] = i, this.data.array[e + 2] = n, this
	},
	setXYZW: function(e, t, i, n, r) {
		return e = e * this.data.stride + this.offset, this.data.array[e + 0] = t, this.data.array[e + 1] = i, this.data.array[e + 2] = n, this.data.array[e + 3] = r, this
	},
	clone: function(e) {
		if (void 0 === e) {
			console.log("THREE.InterleavedBufferAttribute.clone(): Cloning an interlaved buffer attribute will deinterleave buffer data.");
			const e = [];
			for (let t = 0; t < this.count; t++) {
				const i = t * this.data.stride + this.offset;
				for (let t = 0; t < this.itemSize; t++) e.push(this.data.array[i + t])
			}
			return new BufferAttribute(new this.array.constructor(e), this.itemSize, this.normalized)
		}
		return void 0 === e.interleavedBuffers && (e.interleavedBuffers = {}), void 0 === e.interleavedBuffers[this.data.uuid] && (e.interleavedBuffers[this.data.uuid] = this.data.clone(e)), new InterleavedBufferAttribute(e.interleavedBuffers[this.data.uuid], this.itemSize, this.offset, this.normalized)
	},
	toJSON: function(e) {
		if (void 0 === e) {
			console.log("THREE.InterleavedBufferAttribute.toJSON(): Serializing an interlaved buffer attribute will deinterleave buffer data.");
			const e = [];
			for (let t = 0; t < this.count; t++) {
				const i = t * this.data.stride + this.offset;
				for (let t = 0; t < this.itemSize; t++) e.push(this.data.array[i + t])
			}
			return {
				itemSize: this.itemSize,
				type: this.array.constructor.name,
				array: e,
				normalized: this.normalized
			}
		}
		return void 0 === e.interleavedBuffers && (e.interleavedBuffers = {}), void 0 === e.interleavedBuffers[this.data.uuid] && (e.interleavedBuffers[this.data.uuid] = this.data.toJSON(e)), {
			isInterleavedBufferAttribute: !0,
			itemSize: this.itemSize,
			data: this.data.uuid,
			offset: this.offset,
			normalized: this.normalized
		}
	}
});
const _v1$4 = new Vector3,
	_v2$2 = new Vector3;
class LOD extends Object3D {
	constructor() {
		super(), this._currentLevel = 0, this.type = "LOD", Object.defineProperties(this, {
			levels: {
				enumerable: !0,
				value: []
			},
			isLOD: {
				value: !0
			}
		}), this.autoUpdate = !0
	}
	copy(e) {
		super.copy(e, !1);
		const t = e.levels;
		for (let e = 0, i = t.length; e < i; e++) {
			const i = t[e];
			this.addLevel(i.object.clone(), i.distance)
		}
		return this.autoUpdate = e.autoUpdate, this
	}
	addLevel(e, t = 0) {
		t = Math.abs(t);
		const i = this.levels;
		let n;
		for (n = 0; n < i.length && !(t < i[n].distance); n++);
		return i.splice(n, 0, {
			distance: t,
			object: e
		}), this.add(e), this
	}
	getCurrentLevel() {
		return this._currentLevel
	}
	getObjectForDistance(e) {
		const t = this.levels;
		if (t.length > 0) {
			let i, n;
			for (i = 1, n = t.length; i < n && !(e < t[i].distance); i++);
			return t[i - 1].object
		}
		return null
	}
	raycast(e, t) {
		if (this.levels.length > 0) {
			_v1$4.setFromMatrixPosition(this.matrixWorld);
			const i = e.ray.origin.distanceTo(_v1$4);
			this.getObjectForDistance(i).raycast(e, t)
		}
	}
	update(e) {
		const t = this.levels;
		if (t.length > 1) {
			_v1$4.setFromMatrixPosition(e.matrixWorld), _v2$2.setFromMatrixPosition(this.matrixWorld);
			const i = _v1$4.distanceTo(_v2$2) / e.zoom;
			let n, r;
			for (t[0].object.visible = !0, n = 1, r = t.length; n < r && i >= t[n].distance; n++) t[n - 1].object.visible = !1, t[n].object.visible = !0;
			for (this._currentLevel = n - 1; n < r; n++) t[n].object.visible = !1
		}
	}
	toJSON(e) {
		const t = super.toJSON(e);
		!1 === this.autoUpdate && (t.object.autoUpdate = !1), t.object.levels = [];
		const i = this.levels;
		for (let e = 0, n = i.length; e < n; e++) {
			const n = i[e];
			t.object.levels.push({
				object: n.object.uuid,
				distance: n.distance
			})
		}
		return t
	}
}
const _basePosition = new Vector3,
	_skinIndex = new Vector4,
	_skinWeight = new Vector4,
	_vector$7 = new Vector3,
	_matrix$1 = new Matrix4;

function SkinnedMesh(e, t) {
	Mesh.call(this, e, t), this.type = "SkinnedMesh", this.bindMode = "attached", this.bindMatrix = new Matrix4, this.bindMatrixInverse = new Matrix4
}

function Bone() {
	Object3D.call(this), this.type = "Bone"
}
SkinnedMesh.prototype = Object.assign(Object.create(Mesh.prototype), {
	constructor: SkinnedMesh,
	isSkinnedMesh: !0,
	copy: function(e) {
		return Mesh.prototype.copy.call(this, e), this.bindMode = e.bindMode, this.bindMatrix.copy(e.bindMatrix), this.bindMatrixInverse.copy(e.bindMatrixInverse), this.skeleton = e.skeleton, this
	},
	bind: function(e, t) {
		this.skeleton = e, void 0 === t && (this.updateMatrixWorld(!0), this.skeleton.calculateInverses(), t = this.matrixWorld), this.bindMatrix.copy(t), this.bindMatrixInverse.copy(t).invert()
	},
	pose: function() {
		this.skeleton.pose()
	},
	normalizeSkinWeights: function() {
		const e = new Vector4,
			t = this.geometry.attributes.skinWeight;
		for (let i = 0, n = t.count; i < n; i++) {
			e.x = t.getX(i), e.y = t.getY(i), e.z = t.getZ(i), e.w = t.getW(i);
			const n = 1 / e.manhattanLength();
			n !== 1 / 0 ? e.multiplyScalar(n) : e.set(1, 0, 0, 0), t.setXYZW(i, e.x, e.y, e.z, e.w)
		}
	},
	updateMatrixWorld: function(e) {
		Mesh.prototype.updateMatrixWorld.call(this, e), "attached" === this.bindMode ? this.bindMatrixInverse.copy(this.matrixWorld).invert() : "detached" === this.bindMode ? this.bindMatrixInverse.copy(this.bindMatrix).invert() : console.warn("THREE.SkinnedMesh: Unrecognized bindMode: " + this.bindMode)
	},
	boneTransform: function(e, t) {
		const i = this.skeleton,
			n = this.geometry;
		_skinIndex.fromBufferAttribute(n.attributes.skinIndex, e), _skinWeight.fromBufferAttribute(n.attributes.skinWeight, e), _basePosition.fromBufferAttribute(n.attributes.position, e).applyMatrix4(this.bindMatrix), t.set(0, 0, 0);
		for (let e = 0; e < 4; e++) {
			const n = _skinWeight.getComponent(e);
			if (0 !== n) {
				const r = _skinIndex.getComponent(e);
				_matrix$1.multiplyMatrices(i.bones[r].matrixWorld, i.boneInverses[r]), t.addScaledVector(_vector$7.copy(_basePosition).applyMatrix4(_matrix$1), n)
			}
		}
		return t.applyMatrix4(this.bindMatrixInverse)
	}
}), Bone.prototype = Object.assign(Object.create(Object3D.prototype), {
	constructor: Bone,
	isBone: !0
});
const _offsetMatrix = new Matrix4,
	_identityMatrix = new Matrix4;
class Skeleton {
	constructor(e = [], t = []) {
		this.uuid = MathUtils.generateUUID(), this.bones = e.slice(0), this.boneInverses = t, this.boneMatrices = null, this.boneTexture = null, this.boneTextureSize = 0, this.frame = -1, this.init()
	}
	init() {
		const e = this.bones,
			t = this.boneInverses;
		if (this.boneMatrices = new Float32Array(16 * e.length), 0 === t.length) this.calculateInverses();
		else if (e.length !== t.length) {
			console.warn("THREE.Skeleton: Number of inverse bone matrices does not match amount of bones."), this.boneInverses = [];
			for (let e = 0, t = this.bones.length; e < t; e++) this.boneInverses.push(new Matrix4)
		}
	}
	calculateInverses() {
		this.boneInverses.length = 0;
		for (let e = 0, t = this.bones.length; e < t; e++) {
			const t = new Matrix4;
			this.bones[e] && t.copy(this.bones[e].matrixWorld).invert(), this.boneInverses.push(t)
		}
	}
	pose() {
		for (let e = 0, t = this.bones.length; e < t; e++) {
			const t = this.bones[e];
			t && t.matrixWorld.copy(this.boneInverses[e]).invert()
		}
		for (let e = 0, t = this.bones.length; e < t; e++) {
			const t = this.bones[e];
			t && (t.parent && t.parent.isBone ? (t.matrix.copy(t.parent.matrixWorld).invert(), t.matrix.multiply(t.matrixWorld)) : t.matrix.copy(t.matrixWorld), t.matrix.decompose(t.position, t.quaternion, t.scale))
		}
	}
	update() {
		const e = this.bones,
			t = this.boneInverses,
			i = this.boneMatrices,
			n = this.boneTexture;
		for (let n = 0, r = e.length; n < r; n++) {
			const r = e[n] ? e[n].matrixWorld : _identityMatrix;
			_offsetMatrix.multiplyMatrices(r, t[n]), _offsetMatrix.toArray(i, 16 * n)
		}
		null !== n && (n.needsUpdate = !0)
	}
	clone() {
		return new Skeleton(this.bones, this.boneInverses)
	}
	getBoneByName(e) {
		for (let t = 0, i = this.bones.length; t < i; t++) {
			const i = this.bones[t];
			if (i.name === e) return i
		}
	}
	dispose() {
		null !== this.boneTexture && (this.boneTexture.dispose(), this.boneTexture = null)
	}
	fromJSON(e, t) {
		this.uuid = e.uuid;
		for (let i = 0, n = e.bones.length; i < n; i++) {
			const n = e.bones[i];
			let r = t[n];
			void 0 === r && (console.warn("THREE.Skeleton: No bone found with UUID:", n), r = new Bone), this.bones.push(r), this.boneInverses.push((new Matrix4).fromArray(e.boneInverses[i]))
		}
		return this.init(), this
	}
	toJSON() {
		const e = {
			metadata: {
				version: 4.5,
				type: "Skeleton",
				generator: "Skeleton.toJSON"
			},
			bones: [],
			boneInverses: []
		};
		e.uuid = this.uuid;
		const t = this.bones,
			i = this.boneInverses;
		for (let n = 0, r = t.length; n < r; n++) {
			const r = t[n];
			e.bones.push(r.uuid);
			const a = i[n];
			e.boneInverses.push(a.toArray())
		}
		return e
	}
}
const _instanceLocalMatrix = new Matrix4,
	_instanceWorldMatrix = new Matrix4,
	_instanceIntersects = [],
	_mesh$1 = new Mesh;

function InstancedMesh(e, t, i) {
	Mesh.call(this, e, t), this.instanceMatrix = new BufferAttribute(new Float32Array(16 * i), 16), this.instanceColor = null, this.count = i, this.frustumCulled = !1
}
InstancedMesh.prototype = Object.assign(Object.create(Mesh.prototype), {
	constructor: InstancedMesh,
	isInstancedMesh: !0,
	copy: function(e) {
		return Mesh.prototype.copy.call(this, e), this.instanceMatrix.copy(e.instanceMatrix), null !== e.instanceColor && (this.instanceColor = e.instanceColor.clone()), this.count = e.count, this
	},
	getColorAt: function(e, t) {
		t.fromArray(this.instanceColor.array, 3 * e)
	},
	getMatrixAt: function(e, t) {
		t.fromArray(this.instanceMatrix.array, 16 * e)
	},
	raycast: function(e, t) {
		const i = this.matrixWorld,
			n = this.count;
		if (_mesh$1.geometry = this.geometry, _mesh$1.material = this.material, void 0 !== _mesh$1.material)
			for (let r = 0; r < n; r++) {
				this.getMatrixAt(r, _instanceLocalMatrix), _instanceWorldMatrix.multiplyMatrices(i, _instanceLocalMatrix), _mesh$1.matrixWorld = _instanceWorldMatrix, _mesh$1.raycast(e, _instanceIntersects);
				for (let e = 0, i = _instanceIntersects.length; e < i; e++) {
					const i = _instanceIntersects[e];
					i.instanceId = r, i.object = this, t.push(i)
				}
				_instanceIntersects.length = 0
			}
	},
	setColorAt: function(e, t) {
		null === this.instanceColor && (this.instanceColor = new BufferAttribute(new Float32Array(3 * this.count), 3)), t.toArray(this.instanceColor.array, 3 * e)
	},
	setMatrixAt: function(e, t) {
		t.toArray(this.instanceMatrix.array, 16 * e)
	},
	updateMorphTargets: function() {},
	dispose: function() {
		this.dispatchEvent({
			type: "dispose"
		})
	}
});
class LineBasicMaterial extends Material$1 {
	constructor(e) {
		super(), this.type = "LineBasicMaterial", this.color = new Color(16777215), this.linewidth = 1, this.linecap = "round", this.linejoin = "round", this.morphTargets = !1, this.setValues(e)
	}
	copy(e) {
		return super.copy(e), this.color.copy(e.color), this.linewidth = e.linewidth, this.linecap = e.linecap, this.linejoin = e.linejoin, this.morphTargets = e.morphTargets, this
	}
}
LineBasicMaterial.prototype.isLineBasicMaterial = !0;
const _start = new Vector3,
	_end = new Vector3,
	_inverseMatrix$1 = new Matrix4,
	_ray$1 = new Ray,
	_sphere$2 = new Sphere;

function Line(e = new BufferGeometry, t = new LineBasicMaterial) {
	Object3D.call(this), this.type = "Line", this.geometry = e, this.material = t, this.updateMorphTargets()
}
Line.prototype = Object.assign(Object.create(Object3D.prototype), {
	constructor: Line,
	isLine: !0,
	copy: function(e) {
		return Object3D.prototype.copy.call(this, e), this.material = e.material, this.geometry = e.geometry, this
	},
	computeLineDistances: function() {
		const e = this.geometry;
		if (e.isBufferGeometry)
			if (null === e.index) {
				const t = e.attributes.position,
					i = [0];
				for (let e = 1, n = t.count; e < n; e++) _start.fromBufferAttribute(t, e - 1), _end.fromBufferAttribute(t, e), i[e] = i[e - 1], i[e] += _start.distanceTo(_end);
				e.setAttribute("lineDistance", new Float32BufferAttribute(i, 1))
			} else console.warn("THREE.Line.computeLineDistances(): Computation only possible with non-indexed BufferGeometry.");
		else e.isGeometry && console.error("THREE.Line.computeLineDistances() no longer supports THREE.Geometry. Use THREE.BufferGeometry instead.");
		return this
	},
	raycast: function(e, t) {
		const i = this.geometry,
			n = this.matrixWorld,
			r = e.params.Line.threshold;
		if (null === i.boundingSphere && i.computeBoundingSphere(), _sphere$2.copy(i.boundingSphere), _sphere$2.applyMatrix4(n), _sphere$2.radius += r, !1 === e.ray.intersectsSphere(_sphere$2)) return;
		_inverseMatrix$1.copy(n).invert(), _ray$1.copy(e.ray).applyMatrix4(_inverseMatrix$1);
		const a = r / ((this.scale.x + this.scale.y + this.scale.z) / 3),
			s = a * a,
			o = new Vector3,
			l = new Vector3,
			c = new Vector3,
			h = new Vector3,
			u = this.isLineSegments ? 2 : 1;
		if (i.isBufferGeometry) {
			const n = i.index,
				r = i.attributes.position;
			if (null !== n) {
				const i = n.array;
				for (let n = 0, a = i.length - 1; n < a; n += u) {
					const a = i[n],
						u = i[n + 1];
					o.fromBufferAttribute(r, a), l.fromBufferAttribute(r, u);
					if (_ray$1.distanceSqToSegment(o, l, h, c) > s) continue;
					h.applyMatrix4(this.matrixWorld);
					const d = e.ray.origin.distanceTo(h);
					d < e.near || d > e.far || t.push({
						distance: d,
						point: c.clone().applyMatrix4(this.matrixWorld),
						index: n,
						face: null,
						faceIndex: null,
						object: this
					})
				}
			} else
				for (let i = 0, n = r.count - 1; i < n; i += u) {
					o.fromBufferAttribute(r, i), l.fromBufferAttribute(r, i + 1);
					if (_ray$1.distanceSqToSegment(o, l, h, c) > s) continue;
					h.applyMatrix4(this.matrixWorld);
					const n = e.ray.origin.distanceTo(h);
					n < e.near || n > e.far || t.push({
						distance: n,
						point: c.clone().applyMatrix4(this.matrixWorld),
						index: i,
						face: null,
						faceIndex: null,
						object: this
					})
				}
		} else i.isGeometry && console.error("THREE.Line.raycast() no longer supports THREE.Geometry. Use THREE.BufferGeometry instead.")
	},
	updateMorphTargets: function() {
		const e = this.geometry;
		if (e.isBufferGeometry) {
			const t = e.morphAttributes,
				i = Object.keys(t);
			if (i.length > 0) {
				const e = t[i[0]];
				if (void 0 !== e) {
					this.morphTargetInfluences = [], this.morphTargetDictionary = {};
					for (let t = 0, i = e.length; t < i; t++) {
						const i = e[t].name || String(t);
						this.morphTargetInfluences.push(0), this.morphTargetDictionary[i] = t
					}
				}
			}
		} else {
			const t = e.morphTargets;
			void 0 !== t && t.length > 0 && console.error("THREE.Line.updateMorphTargets() does not support THREE.Geometry. Use THREE.BufferGeometry instead.")
		}
	}
});
const _start$1 = new Vector3,
	_end$1 = new Vector3;

function LineSegments(e, t) {
	Line.call(this, e, t), this.type = "LineSegments"
}
LineSegments.prototype = Object.assign(Object.create(Line.prototype), {
	constructor: LineSegments,
	isLineSegments: !0,
	computeLineDistances: function() {
		const e = this.geometry;
		if (e.isBufferGeometry)
			if (null === e.index) {
				const t = e.attributes.position,
					i = [];
				for (let e = 0, n = t.count; e < n; e += 2) _start$1.fromBufferAttribute(t, e), _end$1.fromBufferAttribute(t, e + 1), i[e] = 0 === e ? 0 : i[e - 1], i[e + 1] = i[e] + _start$1.distanceTo(_end$1);
				e.setAttribute("lineDistance", new Float32BufferAttribute(i, 1))
			} else console.warn("THREE.LineSegments.computeLineDistances(): Computation only possible with non-indexed BufferGeometry.");
		else e.isGeometry && console.error("THREE.LineSegments.computeLineDistances() no longer supports THREE.Geometry. Use THREE.BufferGeometry instead.");
		return this
	}
});
class LineLoop extends Line {
	constructor(e, t) {
		super(e, t), this.type = "LineLoop"
	}
}
LineLoop.prototype.isLineLoop = !0;
class PointsMaterial extends Material$1 {
	constructor(e) {
		super(), this.type = "PointsMaterial", this.color = new Color(16777215), this.map = null, this.alphaMap = null, this.size = 1, this.sizeAttenuation = !0, this.morphTargets = !1, this.setValues(e)
	}
	copy(e) {
		return super.copy(e), this.color.copy(e.color), this.map = e.map, this.alphaMap = e.alphaMap, this.size = e.size, this.sizeAttenuation = e.sizeAttenuation, this.morphTargets = e.morphTargets, this
	}
}
PointsMaterial.prototype.isPointsMaterial = !0;
const _inverseMatrix$2 = new Matrix4,
	_ray$2 = new Ray,
	_sphere$3 = new Sphere,
	_position$1 = new Vector3;

function Points(e = new BufferGeometry, t = new PointsMaterial) {
	Object3D.call(this), this.type = "Points", this.geometry = e, this.material = t, this.updateMorphTargets()
}

function testPoint(e, t, i, n, r, a, s) {
	const o = _ray$2.distanceSqToPoint(e);
	if (o < i) {
		const i = new Vector3;
		_ray$2.closestPointToPoint(e, i), i.applyMatrix4(n);
		const l = r.ray.origin.distanceTo(i);
		if (l < r.near || l > r.far) return;
		a.push({
			distance: l,
			distanceToRay: Math.sqrt(o),
			point: i,
			index: t,
			face: null,
			object: s
		})
	}
}
Points.prototype = Object.assign(Object.create(Object3D.prototype), {
	constructor: Points,
	isPoints: !0,
	copy: function(e) {
		return Object3D.prototype.copy.call(this, e), this.material = e.material, this.geometry = e.geometry, this
	},
	raycast: function(e, t) {
		const i = this.geometry,
			n = this.matrixWorld,
			r = e.params.Points.threshold;
		if (null === i.boundingSphere && i.computeBoundingSphere(), _sphere$3.copy(i.boundingSphere), _sphere$3.applyMatrix4(n), _sphere$3.radius += r, !1 === e.ray.intersectsSphere(_sphere$3)) return;
		_inverseMatrix$2.copy(n).invert(), _ray$2.copy(e.ray).applyMatrix4(_inverseMatrix$2);
		const a = r / ((this.scale.x + this.scale.y + this.scale.z) / 3),
			s = a * a;
		if (i.isBufferGeometry) {
			const r = i.index,
				a = i.attributes.position;
			if (null !== r) {
				const i = r.array;
				for (let r = 0, o = i.length; r < o; r++) {
					const o = i[r];
					_position$1.fromBufferAttribute(a, o), testPoint(_position$1, o, s, n, e, t, this)
				}
			} else
				for (let i = 0, r = a.count; i < r; i++) _position$1.fromBufferAttribute(a, i), testPoint(_position$1, i, s, n, e, t, this)
		} else console.error("THREE.Points.raycast() no longer supports THREE.Geometry. Use THREE.BufferGeometry instead.")
	},
	updateMorphTargets: function() {
		const e = this.geometry;
		if (e.isBufferGeometry) {
			const t = e.morphAttributes,
				i = Object.keys(t);
			if (i.length > 0) {
				const e = t[i[0]];
				if (void 0 !== e) {
					this.morphTargetInfluences = [], this.morphTargetDictionary = {};
					for (let t = 0, i = e.length; t < i; t++) {
						const i = e[t].name || String(t);
						this.morphTargetInfluences.push(0), this.morphTargetDictionary[i] = t
					}
				}
			}
		} else {
			const t = e.morphTargets;
			void 0 !== t && t.length > 0 && console.error("THREE.Points.updateMorphTargets() does not support THREE.Geometry. Use THREE.BufferGeometry instead.")
		}
	}
});
class CompressedTexture extends Texture$1 {
	constructor(e, t, i, n, r, a, s, o, l, c, h, u) {
		super(null, a, s, o, l, c, n, r, h, u), this.image = {
			width: t,
			height: i
		}, this.mipmaps = e, this.flipY = !1, this.generateMipmaps = !1
	}
}
CompressedTexture.prototype.isCompressedTexture = !0;
class CanvasTexture extends Texture$1 {
	constructor(e, t, i, n, r, a, s, o, l) {
		super(e, t, i, n, r, a, s, o, l), this.needsUpdate = !0
	}
}

function ParametricGeometry(e, t, i) {
	BufferGeometry.call(this), this.type = "ParametricGeometry", this.parameters = {
		func: e,
		slices: t,
		stacks: i
	};
	const n = [],
		r = [],
		a = [],
		s = [],
		o = 1e-5,
		l = new Vector3,
		c = new Vector3,
		h = new Vector3,
		u = new Vector3,
		d = new Vector3;
	e.length < 3 && console.error("THREE.ParametricGeometry: Function must now modify a Vector3 as third parameter.");
	const p = t + 1;
	for (let n = 0; n <= i; n++) {
		const p = n / i;
		for (let i = 0; i <= t; i++) {
			const n = i / t;
			e(n, p, c), r.push(c.x, c.y, c.z), n - o >= 0 ? (e(n - o, p, h), u.subVectors(c, h)) : (e(n + o, p, h), u.subVectors(h, c)), p - o >= 0 ? (e(n, p - o, h), d.subVectors(c, h)) : (e(n, p + o, h), d.subVectors(h, c)), l.crossVectors(u, d).normalize(), a.push(l.x, l.y, l.z), s.push(n, p)
		}
	}
	for (let e = 0; e < i; e++)
		for (let i = 0; i < t; i++) {
			const t = e * p + i,
				r = e * p + i + 1,
				a = (e + 1) * p + i + 1,
				s = (e + 1) * p + i;
			n.push(t, r, s), n.push(r, a, s)
		}
	this.setIndex(n), this.setAttribute("position", new Float32BufferAttribute(r, 3)), this.setAttribute("normal", new Float32BufferAttribute(a, 3)), this.setAttribute("uv", new Float32BufferAttribute(s, 2))
}
CanvasTexture.prototype.isCanvasTexture = !0, ParametricGeometry.prototype = Object.create(BufferGeometry.prototype), ParametricGeometry.prototype.constructor = ParametricGeometry;
class ShadowMaterial extends Material$1 {
	constructor(e) {
		super(), this.type = "ShadowMaterial", this.color = new Color(0), this.transparent = !0, this.setValues(e)
	}
	copy(e) {
		return super.copy(e), this.color.copy(e.color), this
	}
}
ShadowMaterial.prototype.isShadowMaterial = !0;
class RawShaderMaterial extends ShaderMaterial {
	constructor(e) {
		super(e), this.type = "RawShaderMaterial"
	}
}

function MeshStandardMaterial(e) {
	Material$1.call(this), this.defines = {
		STANDARD: ""
	}, this.type = "MeshStandardMaterial", this.color = new Color(16777215), this.roughness = 1, this.metalness = 0, this.map = null, this.lightMap = null, this.lightMapIntensity = 1, this.aoMap = null, this.aoMapIntensity = 1, this.emissive = new Color(0), this.emissiveIntensity = 1, this.emissiveMap = null, this.bumpMap = null, this.bumpScale = 1, this.normalMap = null, this.normalMapType = 0, this.normalScale = new Vector2(1, 1), this.displacementMap = null, this.displacementScale = 1, this.displacementBias = 0, this.roughnessMap = null, this.metalnessMap = null, this.alphaMap = null, this.envMap = null, this.envMapIntensity = 1, this.refractionRatio = .98, this.wireframe = !1, this.wireframeLinewidth = 1, this.wireframeLinecap = "round", this.wireframeLinejoin = "round", this.skinning = !1, this.morphTargets = !1, this.morphNormals = !1, this.flatShading = !1, this.vertexTangents = !1, this.setValues(e)
}

function MeshPhysicalMaterial(e) {
	MeshStandardMaterial.call(this), this.defines = {
		STANDARD: "",
		PHYSICAL: ""
	}, this.type = "MeshPhysicalMaterial", this.clearcoat = 0, this.clearcoatMap = null, this.clearcoatRoughness = 0, this.clearcoatRoughnessMap = null, this.clearcoatNormalScale = new Vector2(1, 1), this.clearcoatNormalMap = null, this.reflectivity = .5, Object.defineProperty(this, "ior", {
		get: function() {
			return (1 + .4 * this.reflectivity) / (1 - .4 * this.reflectivity)
		},
		set: function(e) {
			this.reflectivity = MathUtils.clamp(2.5 * (e - 1) / (e + 1), 0, 1)
		}
	}), this.sheen = null, this.transmission = 0, this.transmissionMap = null, this.setValues(e)
}
RawShaderMaterial.prototype.isRawShaderMaterial = !0, MeshStandardMaterial.prototype = Object.create(Material$1.prototype), MeshStandardMaterial.prototype.constructor = MeshStandardMaterial, MeshStandardMaterial.prototype.isMeshStandardMaterial = !0, MeshStandardMaterial.prototype.copy = function(e) {
	return Material$1.prototype.copy.call(this, e), this.defines = {
		STANDARD: ""
	}, this.color.copy(e.color), this.roughness = e.roughness, this.metalness = e.metalness, this.map = e.map, this.lightMap = e.lightMap, this.lightMapIntensity = e.lightMapIntensity, this.aoMap = e.aoMap, this.aoMapIntensity = e.aoMapIntensity, this.emissive.copy(e.emissive), this.emissiveMap = e.emissiveMap, this.emissiveIntensity = e.emissiveIntensity, this.bumpMap = e.bumpMap, this.bumpScale = e.bumpScale, this.normalMap = e.normalMap, this.normalMapType = e.normalMapType, this.normalScale.copy(e.normalScale), this.displacementMap = e.displacementMap, this.displacementScale = e.displacementScale, this.displacementBias = e.displacementBias, this.roughnessMap = e.roughnessMap, this.metalnessMap = e.metalnessMap, this.alphaMap = e.alphaMap, this.envMap = e.envMap, this.envMapIntensity = e.envMapIntensity, this.refractionRatio = e.refractionRatio, this.wireframe = e.wireframe, this.wireframeLinewidth = e.wireframeLinewidth, this.wireframeLinecap = e.wireframeLinecap, this.wireframeLinejoin = e.wireframeLinejoin, this.skinning = e.skinning, this.morphTargets = e.morphTargets, this.morphNormals = e.morphNormals, this.flatShading = e.flatShading, this.vertexTangents = e.vertexTangents, this
}, MeshPhysicalMaterial.prototype = Object.create(MeshStandardMaterial.prototype), MeshPhysicalMaterial.prototype.constructor = MeshPhysicalMaterial, MeshPhysicalMaterial.prototype.isMeshPhysicalMaterial = !0, MeshPhysicalMaterial.prototype.copy = function(e) {
	return MeshStandardMaterial.prototype.copy.call(this, e), this.defines = {
		STANDARD: "",
		PHYSICAL: ""
	}, this.clearcoat = e.clearcoat, this.clearcoatMap = e.clearcoatMap, this.clearcoatRoughness = e.clearcoatRoughness, this.clearcoatRoughnessMap = e.clearcoatRoughnessMap, this.clearcoatNormalMap = e.clearcoatNormalMap, this.clearcoatNormalScale.copy(e.clearcoatNormalScale), this.reflectivity = e.reflectivity, e.sheen ? this.sheen = (this.sheen || new Color).copy(e.sheen) : this.sheen = null, this.transmission = e.transmission, this.transmissionMap = e.transmissionMap, this
};
class MeshPhongMaterial extends Material$1 {
	constructor(e) {
		super(), this.type = "MeshPhongMaterial", this.color = new Color(16777215), this.specular = new Color(1118481), this.shininess = 30, this.map = null, this.lightMap = null, this.lightMapIntensity = 1, this.aoMap = null, this.aoMapIntensity = 1, this.emissive = new Color(0), this.emissiveIntensity = 1, this.emissiveMap = null, this.bumpMap = null, this.bumpScale = 1, this.normalMap = null, this.normalMapType = 0, this.normalScale = new Vector2(1, 1), this.displacementMap = null, this.displacementScale = 1, this.displacementBias = 0, this.specularMap = null, this.alphaMap = null, this.envMap = null, this.combine = 0, this.reflectivity = 1, this.refractionRatio = .98, this.wireframe = !1, this.wireframeLinewidth = 1, this.wireframeLinecap = "round", this.wireframeLinejoin = "round", this.skinning = !1, this.morphTargets = !1, this.morphNormals = !1, this.flatShading = !1, this.setValues(e)
	}
	copy(e) {
		return super.copy(e), this.color.copy(e.color), this.specular.copy(e.specular), this.shininess = e.shininess, this.map = e.map, this.lightMap = e.lightMap, this.lightMapIntensity = e.lightMapIntensity, this.aoMap = e.aoMap, this.aoMapIntensity = e.aoMapIntensity, this.emissive.copy(e.emissive), this.emissiveMap = e.emissiveMap, this.emissiveIntensity = e.emissiveIntensity, this.bumpMap = e.bumpMap, this.bumpScale = e.bumpScale, this.normalMap = e.normalMap, this.normalMapType = e.normalMapType, this.normalScale.copy(e.normalScale), this.displacementMap = e.displacementMap, this.displacementScale = e.displacementScale, this.displacementBias = e.displacementBias, this.specularMap = e.specularMap, this.alphaMap = e.alphaMap, this.envMap = e.envMap, this.combine = e.combine, this.reflectivity = e.reflectivity, this.refractionRatio = e.refractionRatio, this.wireframe = e.wireframe, this.wireframeLinewidth = e.wireframeLinewidth, this.wireframeLinecap = e.wireframeLinecap, this.wireframeLinejoin = e.wireframeLinejoin, this.skinning = e.skinning, this.morphTargets = e.morphTargets, this.morphNormals = e.morphNormals, this.flatShading = e.flatShading, this
	}
}
MeshPhongMaterial.prototype.isMeshPhongMaterial = !0;
const AnimationUtils = {
	arraySlice: function(e, t, i) {
		return AnimationUtils.isTypedArray(e) ? new e.constructor(e.subarray(t, void 0 !== i ? i : e.length)) : e.slice(t, i)
	},
	convertArray: function(e, t, i) {
		return !e || !i && e.constructor === t ? e : "number" == typeof t.BYTES_PER_ELEMENT ? new t(e) : Array.prototype.slice.call(e)
	},
	isTypedArray: function(e) {
		return ArrayBuffer.isView(e) && !(e instanceof DataView)
	},
	getKeyframeOrder: function(e) {
		const t = e.length,
			i = new Array(t);
		for (let e = 0; e !== t; ++e) i[e] = e;
		return i.sort((function(t, i) {
			return e[t] - e[i]
		})), i
	},
	sortedArray: function(e, t, i) {
		const n = e.length,
			r = new e.constructor(n);
		for (let a = 0, s = 0; s !== n; ++a) {
			const n = i[a] * t;
			for (let i = 0; i !== t; ++i) r[s++] = e[n + i]
		}
		return r
	},
	flattenJSON: function(e, t, i, n) {
		let r = 1,
			a = e[0];
		for (; void 0 !== a && void 0 === a[n];) a = e[r++];
		if (void 0 === a) return;
		let s = a[n];
		if (void 0 !== s)
			if (Array.isArray(s))
				do {
					s = a[n], void 0 !== s && (t.push(a.time), i.push.apply(i, s)), a = e[r++]
				} while (void 0 !== a);
			else if (void 0 !== s.toArray)
			do {
				s = a[n], void 0 !== s && (t.push(a.time), s.toArray(i, i.length)), a = e[r++]
			} while (void 0 !== a);
		else
			do {
				s = a[n], void 0 !== s && (t.push(a.time), i.push(s)), a = e[r++]
			} while (void 0 !== a)
	},
	subclip: function(e, t, i, n, r = 30) {
		const a = e.clone();
		a.name = t;
		const s = [];
		for (let e = 0; e < a.tracks.length; ++e) {
			const t = a.tracks[e],
				o = t.getValueSize(),
				l = [],
				c = [];
			for (let e = 0; e < t.times.length; ++e) {
				const a = t.times[e] * r;
				if (!(a < i || a >= n)) {
					l.push(t.times[e]);
					for (let i = 0; i < o; ++i) c.push(t.values[e * o + i])
				}
			}
			0 !== l.length && (t.times = AnimationUtils.convertArray(l, t.times.constructor), t.values = AnimationUtils.convertArray(c, t.values.constructor), s.push(t))
		}
		a.tracks = s;
		let o = 1 / 0;
		for (let e = 0; e < a.tracks.length; ++e) o > a.tracks[e].times[0] && (o = a.tracks[e].times[0]);
		for (let e = 0; e < a.tracks.length; ++e) a.tracks[e].shift(-1 * o);
		return a.resetDuration(), a
	},
	makeClipAdditive: function(e, t = 0, i = e, n = 30) {
		n <= 0 && (n = 30);
		const r = i.tracks.length,
			a = t / n;
		for (let t = 0; t < r; ++t) {
			const n = i.tracks[t],
				r = n.ValueTypeName;
			if ("bool" === r || "string" === r) continue;
			const s = e.tracks.find((function(e) {
				return e.name === n.name && e.ValueTypeName === r
			}));
			if (void 0 === s) continue;
			let o = 0;
			const l = n.getValueSize();
			n.createInterpolant.isInterpolantFactoryMethodGLTFCubicSpline && (o = l / 3);
			let c = 0;
			const h = s.getValueSize();
			s.createInterpolant.isInterpolantFactoryMethodGLTFCubicSpline && (c = h / 3);
			const u = n.times.length - 1;
			let d;
			if (a <= n.times[0]) {
				const e = o,
					t = l - o;
				d = AnimationUtils.arraySlice(n.values, e, t)
			} else if (a >= n.times[u]) {
				const e = u * l + o,
					t = e + l - o;
				d = AnimationUtils.arraySlice(n.values, e, t)
			} else {
				const e = n.createInterpolant(),
					t = o,
					i = l - o;
				e.evaluate(a), d = AnimationUtils.arraySlice(e.resultBuffer, t, i)
			}
			if ("quaternion" === r) {
				(new Quaternion).fromArray(d).normalize().conjugate().toArray(d)
			}
			const p = s.times.length;
			for (let e = 0; e < p; ++e) {
				const t = e * h + c;
				if ("quaternion" === r) Quaternion.multiplyQuaternionsFlat(s.values, t, d, 0, s.values, t);
				else {
					const e = h - 2 * c;
					for (let i = 0; i < e; ++i) s.values[t + i] -= d[i]
				}
			}
		}
		return e.blendMode = 2501, e
	}
};

function Interpolant(e, t, i, n) {
	this.parameterPositions = e, this._cachedIndex = 0, this.resultBuffer = void 0 !== n ? n : new t.constructor(i), this.sampleValues = t, this.valueSize = i
}

function CubicInterpolant(e, t, i, n) {
	Interpolant.call(this, e, t, i, n), this._weightPrev = -0, this._offsetPrev = -0, this._weightNext = -0, this._offsetNext = -0
}

function LinearInterpolant(e, t, i, n) {
	Interpolant.call(this, e, t, i, n)
}

function DiscreteInterpolant(e, t, i, n) {
	Interpolant.call(this, e, t, i, n)
}
Object.assign(Interpolant.prototype, {
	evaluate: function(e) {
		const t = this.parameterPositions;
		let i = this._cachedIndex,
			n = t[i],
			r = t[i - 1];
		e: {
			t: {
				let a;i: {
					n: if (!(e < n)) {
						for (let a = i + 2;;) {
							if (void 0 === n) {
								if (e < r) break n;
								return i = t.length, this._cachedIndex = i, this.afterEnd_(i - 1, e, r)
							}
							if (i === a) break;
							if (r = n, n = t[++i], e < n) break t
						}
						a = t.length;
						break i
					}if (e >= r) break e; {
						const s = t[1];
						e < s && (i = 2, r = s);
						for (let a = i - 2;;) {
							if (void 0 === r) return this._cachedIndex = 0, this.beforeStart_(0, e, n);
							if (i === a) break;
							if (n = r, r = t[--i - 1], e >= r) break t
						}
						a = i, i = 0
					}
				}
				for (; i < a;) {
					const n = i + a >>> 1;
					e < t[n] ? a = n : i = n + 1
				}
				if (n = t[i], r = t[i - 1], void 0 === r) return this._cachedIndex = 0,
				this.beforeStart_(0, e, n);
				if (void 0 === n) return i = t.length,
				this._cachedIndex = i,
				this.afterEnd_(i - 1, r, e)
			}
			this._cachedIndex = i,
			this.intervalChanged_(i, r, n)
		}
		return this.interpolate_(i, r, e, n)
	},
	settings: null,
	DefaultSettings_: {},
	getSettings_: function() {
		return this.settings || this.DefaultSettings_
	},
	copySampleValue_: function(e) {
		const t = this.resultBuffer,
			i = this.sampleValues,
			n = this.valueSize,
			r = e * n;
		for (let e = 0; e !== n; ++e) t[e] = i[r + e];
		return t
	},
	interpolate_: function() {
		throw new Error("call to abstract method")
	},
	intervalChanged_: function() {}
}), Object.assign(Interpolant.prototype, {
	beforeStart_: Interpolant.prototype.copySampleValue_,
	afterEnd_: Interpolant.prototype.copySampleValue_
}), CubicInterpolant.prototype = Object.assign(Object.create(Interpolant.prototype), {
	constructor: CubicInterpolant,
	DefaultSettings_: {
		endingStart: 2400,
		endingEnd: 2400
	},
	intervalChanged_: function(e, t, i) {
		const n = this.parameterPositions;
		let r = e - 2,
			a = e + 1,
			s = n[r],
			o = n[a];
		if (void 0 === s) switch (this.getSettings_().endingStart) {
			case 2401:
				r = e, s = 2 * t - i;
				break;
			case 2402:
				r = n.length - 2, s = t + n[r] - n[r + 1];
				break;
			default:
				r = e, s = i
		}
		if (void 0 === o) switch (this.getSettings_().endingEnd) {
			case 2401:
				a = e, o = 2 * i - t;
				break;
			case 2402:
				a = 1, o = i + n[1] - n[0];
				break;
			default:
				a = e - 1, o = t
		}
		const l = .5 * (i - t),
			c = this.valueSize;
		this._weightPrev = l / (t - s), this._weightNext = l / (o - i), this._offsetPrev = r * c, this._offsetNext = a * c
	},
	interpolate_: function(e, t, i, n) {
		const r = this.resultBuffer,
			a = this.sampleValues,
			s = this.valueSize,
			o = e * s,
			l = o - s,
			c = this._offsetPrev,
			h = this._offsetNext,
			u = this._weightPrev,
			d = this._weightNext,
			p = (i - t) / (n - t),
			m = p * p,
			A = m * p,
			g = -u * A + 2 * u * m - u * p,
			f = (1 + u) * A + (-1.5 - 2 * u) * m + (-.5 + u) * p + 1,
			v = (-1 - d) * A + (1.5 + d) * m + .5 * p,
			y = d * A - d * m;
		for (let e = 0; e !== s; ++e) r[e] = g * a[c + e] + f * a[l + e] + v * a[o + e] + y * a[h + e];
		return r
	}
}), LinearInterpolant.prototype = Object.assign(Object.create(Interpolant.prototype), {
	constructor: LinearInterpolant,
	interpolate_: function(e, t, i, n) {
		const r = this.resultBuffer,
			a = this.sampleValues,
			s = this.valueSize,
			o = e * s,
			l = o - s,
			c = (i - t) / (n - t),
			h = 1 - c;
		for (let e = 0; e !== s; ++e) r[e] = a[l + e] * h + a[o + e] * c;
		return r
	}
}), DiscreteInterpolant.prototype = Object.assign(Object.create(Interpolant.prototype), {
	constructor: DiscreteInterpolant,
	interpolate_: function(e) {
		return this.copySampleValue_(e - 1)
	}
});
class KeyframeTrack {
	constructor(e, t, i, n) {
		if (void 0 === e) throw new Error("THREE.KeyframeTrack: track name is undefined");
		if (void 0 === t || 0 === t.length) throw new Error("THREE.KeyframeTrack: no keyframes in track named " + e);
		this.name = e, this.times = AnimationUtils.convertArray(t, this.TimeBufferType), this.values = AnimationUtils.convertArray(i, this.ValueBufferType), this.setInterpolation(n || this.DefaultInterpolation)
	}
	static toJSON(e) {
		const t = e.constructor;
		let i;
		if (t.toJSON !== this.toJSON) i = t.toJSON(e);
		else {
			i = {
				name: e.name,
				times: AnimationUtils.convertArray(e.times, Array),
				values: AnimationUtils.convertArray(e.values, Array)
			};
			const t = e.getInterpolation();
			t !== e.DefaultInterpolation && (i.interpolation = t)
		}
		return i.type = e.ValueTypeName, i
	}
	InterpolantFactoryMethodDiscrete(e) {
		return new DiscreteInterpolant(this.times, this.values, this.getValueSize(), e)
	}
	InterpolantFactoryMethodLinear(e) {
		return new LinearInterpolant(this.times, this.values, this.getValueSize(), e)
	}
	InterpolantFactoryMethodSmooth(e) {
		return new CubicInterpolant(this.times, this.values, this.getValueSize(), e)
	}
	setInterpolation(e) {
		let t;
		switch (e) {
			case 2300:
				t = this.InterpolantFactoryMethodDiscrete;
				break;
			case 2301:
				t = this.InterpolantFactoryMethodLinear;
				break;
			case 2302:
				t = this.InterpolantFactoryMethodSmooth
		}
		if (void 0 === t) {
			const t = "unsupported interpolation for " + this.ValueTypeName + " keyframe track named " + this.name;
			if (void 0 === this.createInterpolant) {
				if (e === this.DefaultInterpolation) throw new Error(t);
				this.setInterpolation(this.DefaultInterpolation)
			}
			return console.warn("THREE.KeyframeTrack:", t), this
		}
		return this.createInterpolant = t, this
	}
	getInterpolation() {
		switch (this.createInterpolant) {
			case this.InterpolantFactoryMethodDiscrete:
				return 2300;
			case this.InterpolantFactoryMethodLinear:
				return 2301;
			case this.InterpolantFactoryMethodSmooth:
				return 2302
		}
	}
	getValueSize() {
		return this.values.length / this.times.length
	}
	shift(e) {
		if (0 !== e) {
			const t = this.times;
			for (let i = 0, n = t.length; i !== n; ++i) t[i] += e
		}
		return this
	}
	scale(e) {
		if (1 !== e) {
			const t = this.times;
			for (let i = 0, n = t.length; i !== n; ++i) t[i] *= e
		}
		return this
	}
	trim(e, t) {
		const i = this.times,
			n = i.length;
		let r = 0,
			a = n - 1;
		for (; r !== n && i[r] < e;) ++r;
		for (; - 1 !== a && i[a] > t;) --a;
		if (++a, 0 !== r || a !== n) {
			r >= a && (a = Math.max(a, 1), r = a - 1);
			const e = this.getValueSize();
			this.times = AnimationUtils.arraySlice(i, r, a), this.values = AnimationUtils.arraySlice(this.values, r * e, a * e)
		}
		return this
	}
	validate() {
		let e = !0;
		const t = this.getValueSize();
		t - Math.floor(t) != 0 && (console.error("THREE.KeyframeTrack: Invalid value size in track.", this), e = !1);
		const i = this.times,
			n = this.values,
			r = i.length;
		0 === r && (console.error("THREE.KeyframeTrack: Track is empty.", this), e = !1);
		let a = null;
		for (let t = 0; t !== r; t++) {
			const n = i[t];
			if ("number" == typeof n && isNaN(n)) {
				console.error("THREE.KeyframeTrack: Time is not a valid number.", this, t, n), e = !1;
				break
			}
			if (null !== a && a > n) {
				console.error("THREE.KeyframeTrack: Out of order keys.", this, t, n, a), e = !1;
				break
			}
			a = n
		}
		if (void 0 !== n && AnimationUtils.isTypedArray(n))
			for (let t = 0, i = n.length; t !== i; ++t) {
				const i = n[t];
				if (isNaN(i)) {
					console.error("THREE.KeyframeTrack: Value is not a valid number.", this, t, i), e = !1;
					break
				}
			}
		return e
	}
	optimize() {
		const e = AnimationUtils.arraySlice(this.times),
			t = AnimationUtils.arraySlice(this.values),
			i = this.getValueSize(),
			n = 2302 === this.getInterpolation(),
			r = e.length - 1;
		let a = 1;
		for (let s = 1; s < r; ++s) {
			let r = !1;
			const o = e[s];
			if (o !== e[s + 1] && (1 !== s || o !== e[0]))
				if (n) r = !0;
				else {
					const e = s * i,
						n = e - i,
						a = e + i;
					for (let s = 0; s !== i; ++s) {
						const i = t[e + s];
						if (i !== t[n + s] || i !== t[a + s]) {
							r = !0;
							break
						}
					}
				} if (r) {
				if (s !== a) {
					e[a] = e[s];
					const n = s * i,
						r = a * i;
					for (let e = 0; e !== i; ++e) t[r + e] = t[n + e]
				}++a
			}
		}
		if (r > 0) {
			e[a] = e[r];
			for (let e = r * i, n = a * i, s = 0; s !== i; ++s) t[n + s] = t[e + s];
			++a
		}
		return a !== e.length ? (this.times = AnimationUtils.arraySlice(e, 0, a), this.values = AnimationUtils.arraySlice(t, 0, a * i)) : (this.times = e, this.values = t), this
	}
	clone() {
		const e = AnimationUtils.arraySlice(this.times, 0),
			t = AnimationUtils.arraySlice(this.values, 0),
			i = new(0, this.constructor)(this.name, e, t);
		return i.createInterpolant = this.createInterpolant, i
	}
}
KeyframeTrack.prototype.TimeBufferType = Float32Array, KeyframeTrack.prototype.ValueBufferType = Float32Array, KeyframeTrack.prototype.DefaultInterpolation = 2301;
class BooleanKeyframeTrack extends KeyframeTrack {}
BooleanKeyframeTrack.prototype.ValueTypeName = "bool", BooleanKeyframeTrack.prototype.ValueBufferType = Array, BooleanKeyframeTrack.prototype.DefaultInterpolation = 2300, BooleanKeyframeTrack.prototype.InterpolantFactoryMethodLinear = void 0, BooleanKeyframeTrack.prototype.InterpolantFactoryMethodSmooth = void 0;
class ColorKeyframeTrack extends KeyframeTrack {}
ColorKeyframeTrack.prototype.ValueTypeName = "color";
class NumberKeyframeTrack extends KeyframeTrack {}

function QuaternionLinearInterpolant(e, t, i, n) {
	Interpolant.call(this, e, t, i, n)
}
NumberKeyframeTrack.prototype.ValueTypeName = "number", QuaternionLinearInterpolant.prototype = Object.assign(Object.create(Interpolant.prototype), {
	constructor: QuaternionLinearInterpolant,
	interpolate_: function(e, t, i, n) {
		const r = this.resultBuffer,
			a = this.sampleValues,
			s = this.valueSize,
			o = (i - t) / (n - t);
		let l = e * s;
		for (let e = l + s; l !== e; l += 4) Quaternion.slerpFlat(r, 0, a, l - s, a, l, o);
		return r
	}
});
class QuaternionKeyframeTrack extends KeyframeTrack {
	InterpolantFactoryMethodLinear(e) {
		return new QuaternionLinearInterpolant(this.times, this.values, this.getValueSize(), e)
	}
}
QuaternionKeyframeTrack.prototype.ValueTypeName = "quaternion", QuaternionKeyframeTrack.prototype.DefaultInterpolation = 2301, QuaternionKeyframeTrack.prototype.InterpolantFactoryMethodSmooth = void 0;
class StringKeyframeTrack extends KeyframeTrack {}
StringKeyframeTrack.prototype.ValueTypeName = "string", StringKeyframeTrack.prototype.ValueBufferType = Array, StringKeyframeTrack.prototype.DefaultInterpolation = 2300, StringKeyframeTrack.prototype.InterpolantFactoryMethodLinear = void 0, StringKeyframeTrack.prototype.InterpolantFactoryMethodSmooth = void 0;
class VectorKeyframeTrack extends KeyframeTrack {}
VectorKeyframeTrack.prototype.ValueTypeName = "vector";
class AnimationClip {
	constructor(e, t = -1, i, n = 2500) {
		this.name = e, this.tracks = i, this.duration = t, this.blendMode = n, this.uuid = MathUtils.generateUUID(), this.duration < 0 && this.resetDuration()
	}
	static parse(e) {
		const t = [],
			i = e.tracks,
			n = 1 / (e.fps || 1);
		for (let e = 0, r = i.length; e !== r; ++e) t.push(parseKeyframeTrack(i[e]).scale(n));
		const r = new this(e.name, e.duration, t, e.blendMode);
		return r.uuid = e.uuid, r
	}
	static toJSON(e) {
		const t = [],
			i = e.tracks,
			n = {
				name: e.name,
				duration: e.duration,
				tracks: t,
				uuid: e.uuid,
				blendMode: e.blendMode
			};
		for (let e = 0, n = i.length; e !== n; ++e) t.push(KeyframeTrack.toJSON(i[e]));
		return n
	}
	static CreateFromMorphTargetSequence(e, t, i, n) {
		const r = t.length,
			a = [];
		for (let e = 0; e < r; e++) {
			let s = [],
				o = [];
			s.push((e + r - 1) % r, e, (e + 1) % r), o.push(0, 1, 0);
			const l = AnimationUtils.getKeyframeOrder(s);
			s = AnimationUtils.sortedArray(s, 1, l), o = AnimationUtils.sortedArray(o, 1, l), n || 0 !== s[0] || (s.push(r), o.push(o[0])), a.push(new NumberKeyframeTrack(".morphTargetInfluences[" + t[e].name + "]", s, o).scale(1 / i))
		}
		return new this(e, -1, a)
	}
	static findByName(e, t) {
		let i = e;
		if (!Array.isArray(e)) {
			const t = e;
			i = t.geometry && t.geometry.animations || t.animations
		}
		for (let e = 0; e < i.length; e++)
			if (i[e].name === t) return i[e];
		return null
	}
	static CreateClipsFromMorphTargetSequences(e, t, i) {
		const n = {},
			r = /^([\w-]*?)([\d]+)$/;
		for (let t = 0, i = e.length; t < i; t++) {
			const i = e[t],
				a = i.name.match(r);
			if (a && a.length > 1) {
				const e = a[1];
				let t = n[e];
				t || (n[e] = t = []), t.push(i)
			}
		}
		const a = [];
		for (const e in n) a.push(this.CreateFromMorphTargetSequence(e, n[e], t, i));
		return a
	}
	static parseAnimation(e, t) {
		if (!e) return console.error("THREE.AnimationClip: No animation in JSONLoader data."), null;
		const i = function(e, t, i, n, r) {
				if (0 !== i.length) {
					const a = [],
						s = [];
					AnimationUtils.flattenJSON(i, a, s, n), 0 !== a.length && r.push(new e(t, a, s))
				}
			},
			n = [],
			r = e.name || "default",
			a = e.fps || 30,
			s = e.blendMode;
		let o = e.length || -1;
		const l = e.hierarchy || [];
		for (let e = 0; e < l.length; e++) {
			const r = l[e].keys;
			if (r && 0 !== r.length)
				if (r[0].morphTargets) {
					const e = {};
					let t;
					for (t = 0; t < r.length; t++)
						if (r[t].morphTargets)
							for (let i = 0; i < r[t].morphTargets.length; i++) e[r[t].morphTargets[i]] = -1;
					for (const i in e) {
						const e = [],
							a = [];
						for (let n = 0; n !== r[t].morphTargets.length; ++n) {
							const n = r[t];
							e.push(n.time), a.push(n.morphTarget === i ? 1 : 0)
						}
						n.push(new NumberKeyframeTrack(".morphTargetInfluence[" + i + "]", e, a))
					}
					o = e.length * (a || 1)
				} else {
					const a = ".bones[" + t[e].name + "]";
					i(VectorKeyframeTrack, a + ".position", r, "pos", n), i(QuaternionKeyframeTrack, a + ".quaternion", r, "rot", n), i(VectorKeyframeTrack, a + ".scale", r, "scl", n)
				}
		}
		if (0 === n.length) return null;
		return new this(r, o, n, s)
	}
	resetDuration() {
		let e = 0;
		for (let t = 0, i = this.tracks.length; t !== i; ++t) {
			const i = this.tracks[t];
			e = Math.max(e, i.times[i.times.length - 1])
		}
		return this.duration = e, this
	}
	trim() {
		for (let e = 0; e < this.tracks.length; e++) this.tracks[e].trim(0, this.duration);
		return this
	}
	validate() {
		let e = !0;
		for (let t = 0; t < this.tracks.length; t++) e = e && this.tracks[t].validate();
		return e
	}
	optimize() {
		for (let e = 0; e < this.tracks.length; e++) this.tracks[e].optimize();
		return this
	}
	clone() {
		const e = [];
		for (let t = 0; t < this.tracks.length; t++) e.push(this.tracks[t].clone());
		return new this.constructor(this.name, this.duration, e, this.blendMode)
	}
	toJSON() {
		return this.constructor.toJSON(this)
	}
}

function getTrackTypeForValueTypeName(e) {
	switch (e.toLowerCase()) {
		case "scalar":
		case "double":
		case "float":
		case "number":
		case "integer":
			return NumberKeyframeTrack;
		case "vector":
		case "vector2":
		case "vector3":
		case "vector4":
			return VectorKeyframeTrack;
		case "color":
			return ColorKeyframeTrack;
		case "quaternion":
			return QuaternionKeyframeTrack;
		case "bool":
		case "boolean":
			return BooleanKeyframeTrack;
		case "string":
			return StringKeyframeTrack
	}
	throw new Error("THREE.KeyframeTrack: Unsupported typeName: " + e)
}

function parseKeyframeTrack(e) {
	if (void 0 === e.type) throw new Error("THREE.KeyframeTrack: track type undefined, can not parse");
	const t = getTrackTypeForValueTypeName(e.type);
	if (void 0 === e.times) {
		const t = [],
			i = [];
		AnimationUtils.flattenJSON(e.keys, t, i, "value"), e.times = t, e.values = i
	}
	return void 0 !== t.parse ? t.parse(e) : new t(e.name, e.times, e.values, e.interpolation)
}
const Cache = {
	enabled: !1,
	files: {},
	add: function(e, t) {
		!1 !== this.enabled && (this.files[e] = t)
	},
	get: function(e) {
		if (!1 !== this.enabled) return this.files[e]
	},
	remove: function(e) {
		delete this.files[e]
	},
	clear: function() {
		this.files = {}
	}
};

function LoadingManager(e, t, i) {
	const n = this;
	let r = !1,
		a = 0,
		s = 0,
		o = void 0;
	const l = [];
	this.onStart = void 0, this.onLoad = e, this.onProgress = t, this.onError = i, this.itemStart = function(e) {
		s++, !1 === r && void 0 !== n.onStart && n.onStart(e, a, s), r = !0
	}, this.itemEnd = function(e) {
		a++, void 0 !== n.onProgress && n.onProgress(e, a, s), a === s && (r = !1, void 0 !== n.onLoad && n.onLoad())
	}, this.itemError = function(e) {
		void 0 !== n.onError && n.onError(e)
	}, this.resolveURL = function(e) {
		return o ? o(e) : e
	}, this.setURLModifier = function(e) {
		return o = e, this
	}, this.addHandler = function(e, t) {
		return l.push(e, t), this
	}, this.removeHandler = function(e) {
		const t = l.indexOf(e);
		return -1 !== t && l.splice(t, 2), this
	}, this.getHandler = function(e) {
		for (let t = 0, i = l.length; t < i; t += 2) {
			const i = l[t],
				n = l[t + 1];
			if (i.global && (i.lastIndex = 0), i.test(e)) return n
		}
		return null
	}
}
const DefaultLoadingManager = new LoadingManager;

function Loader(e) {
	this.manager = void 0 !== e ? e : DefaultLoadingManager, this.crossOrigin = "anonymous", this.withCredentials = !1, this.path = "", this.resourcePath = "", this.requestHeader = {}
}
Object.assign(Loader.prototype, {
	load: function() {},
	loadAsync: function(e, t) {
		const i = this;
		return new Promise((function(n, r) {
			i.load(e, n, t, r)
		}))
	},
	parse: function() {},
	setCrossOrigin: function(e) {
		return this.crossOrigin = e, this
	},
	setWithCredentials: function(e) {
		return this.withCredentials = e, this
	},
	setPath: function(e) {
		return this.path = e, this
	},
	setResourcePath: function(e) {
		return this.resourcePath = e, this
	},
	setRequestHeader: function(e) {
		return this.requestHeader = e, this
	}
});
const loading = {};

function FileLoader(e) {
	Loader.call(this, e)
}

function CompressedTextureLoader(e) {
	Loader.call(this, e)
}
FileLoader.prototype = Object.assign(Object.create(Loader.prototype), {
	constructor: FileLoader,
	load: function(e, t, i, n) {
		void 0 === e && (e = ""), void 0 !== this.path && (e = this.path + e), e = this.manager.resolveURL(e);
		const r = this,
			a = Cache.get(e);
		if (void 0 !== a) return r.manager.itemStart(e), setTimeout((function() {
			t && t(a), r.manager.itemEnd(e)
		}), 0), a;
		if (void 0 !== loading[e]) return void loading[e].push({
			onLoad: t,
			onProgress: i,
			onError: n
		});
		const s = e.match(/^data:(.*?)(;base64)?,(.*)$/);
		let o;
		if (s) {
			const i = s[1],
				a = !!s[2];
			let o = s[3];
			o = decodeURIComponent(o), a && (o = atob(o));
			try {
				let n;
				const a = (this.responseType || "").toLowerCase();
				switch (a) {
					case "arraybuffer":
					case "blob":
						const e = new Uint8Array(o.length);
						for (let t = 0; t < o.length; t++) e[t] = o.charCodeAt(t);
						n = "blob" === a ? new Blob([e.buffer], {
							type: i
						}) : e.buffer;
						break;
					case "document":
						const t = new DOMParser;
						n = t.parseFromString(o, i);
						break;
					case "json":
						n = JSON.parse(o);
						break;
					default:
						n = o
				}
				setTimeout((function() {
					t && t(n), r.manager.itemEnd(e)
				}), 0)
			} catch (t) {
				setTimeout((function() {
					n && n(t), r.manager.itemError(e), r.manager.itemEnd(e)
				}), 0)
			}
		} else {
			loading[e] = [], loading[e].push({
				onLoad: t,
				onProgress: i,
				onError: n
			}), o = new XMLHttpRequest, o.open("GET", e, !0), o.addEventListener("load", (function(t) {
				const i = this.response,
					n = loading[e];
				if (delete loading[e], 200 === this.status || 0 === this.status) {
					0 === this.status && console.warn("THREE.FileLoader: HTTP Status 0 received."), Cache.add(e, i);
					for (let e = 0, t = n.length; e < t; e++) {
						const t = n[e];
						t.onLoad && t.onLoad(i)
					}
					r.manager.itemEnd(e)
				} else {
					for (let e = 0, i = n.length; e < i; e++) {
						const i = n[e];
						i.onError && i.onError(t)
					}
					r.manager.itemError(e), r.manager.itemEnd(e)
				}
			}), !1), o.addEventListener("progress", (function(t) {
				const i = loading[e];
				for (let e = 0, n = i.length; e < n; e++) {
					const n = i[e];
					n.onProgress && n.onProgress(t)
				}
			}), !1), o.addEventListener("error", (function(t) {
				const i = loading[e];
				delete loading[e];
				for (let e = 0, n = i.length; e < n; e++) {
					const n = i[e];
					n.onError && n.onError(t)
				}
				r.manager.itemError(e), r.manager.itemEnd(e)
			}), !1), o.addEventListener("abort", (function(t) {
				const i = loading[e];
				delete loading[e];
				for (let e = 0, n = i.length; e < n; e++) {
					const n = i[e];
					n.onError && n.onError(t)
				}
				r.manager.itemError(e), r.manager.itemEnd(e)
			}), !1), void 0 !== this.responseType && (o.responseType = this.responseType), void 0 !== this.withCredentials && (o.withCredentials = this.withCredentials), o.overrideMimeType && o.overrideMimeType(void 0 !== this.mimeType ? this.mimeType : "text/plain");
			for (const e in this.requestHeader) o.setRequestHeader(e, this.requestHeader[e]);
			o.send(null)
		}
		return r.manager.itemStart(e), o
	},
	setResponseType: function(e) {
		return this.responseType = e, this
	},
	setMimeType: function(e) {
		return this.mimeType = e, this
	}
}), CompressedTextureLoader.prototype = Object.assign(Object.create(Loader.prototype), {
	constructor: CompressedTextureLoader,
	load: function(e, t, i, n) {
		const r = this,
			a = [],
			s = new CompressedTexture,
			o = new FileLoader(this.manager);
		o.setPath(this.path), o.setResponseType("arraybuffer"), o.setRequestHeader(this.requestHeader), o.setWithCredentials(r.withCredentials);
		let l = 0;

		function c(c) {
			o.load(e[c], (function(e) {
				const i = r.parse(e, !0);
				a[c] = {
					width: i.width,
					height: i.height,
					format: i.format,
					mipmaps: i.mipmaps
				}, l += 1, 6 === l && (1 === i.mipmapCount && (s.minFilter = 1006), s.image = a, s.format = i.format, s.needsUpdate = !0, t && t(s))
			}), i, n)
		}
		if (Array.isArray(e))
			for (let t = 0, i = e.length; t < i; ++t) c(t);
		else o.load(e, (function(e) {
			const i = r.parse(e, !0);
			if (i.isCubemap) {
				const e = i.mipmaps.length / i.mipmapCount;
				for (let t = 0; t < e; t++) {
					a[t] = {
						mipmaps: []
					};
					for (let e = 0; e < i.mipmapCount; e++) a[t].mipmaps.push(i.mipmaps[t * i.mipmapCount + e]), a[t].format = i.format, a[t].width = i.width, a[t].height = i.height
				}
				s.image = a
			} else s.image.width = i.width, s.image.height = i.height, s.mipmaps = i.mipmaps;
			1 === i.mipmapCount && (s.minFilter = 1006), s.format = i.format, s.needsUpdate = !0, t && t(s)
		}), i, n);
		return s
	}
});
class ImageLoader extends Loader {
	constructor(e) {
		super(e)
	}
	load(e, t, i, n) {
		void 0 !== this.path && (e = this.path + e), e = this.manager.resolveURL(e);
		const r = this,
			a = Cache.get(e);
		if (void 0 !== a) return r.manager.itemStart(e), setTimeout((function() {
			t && t(a), r.manager.itemEnd(e)
		}), 0), a;
		const s = document.createElementNS("http://www.w3.org/1999/xhtml", "img");

		function o() {
			s.removeEventListener("load", o, !1), s.removeEventListener("error", l, !1), Cache.add(e, this), t && t(this), r.manager.itemEnd(e)
		}

		function l(t) {
			s.removeEventListener("load", o, !1), s.removeEventListener("error", l, !1), n && n(t), r.manager.itemError(e), r.manager.itemEnd(e)
		}
		return s.addEventListener("load", o, !1), s.addEventListener("error", l, !1), "data:" !== e.substr(0, 5) && void 0 !== this.crossOrigin && (s.crossOrigin = this.crossOrigin), r.manager.itemStart(e), s.src = e, s
	}
}
class CubeTextureLoader extends Loader {
	constructor(e) {
		super(e)
	}
	load(e, t, i, n) {
		const r = new CubeTexture,
			a = new ImageLoader(this.manager);
		a.setCrossOrigin(this.crossOrigin), a.setPath(this.path);
		let s = 0;

		function o(i) {
			a.load(e[i], (function(e) {
				r.images[i] = e, s++, 6 === s && (r.needsUpdate = !0, t && t(r))
			}), void 0, n)
		}
		for (let t = 0; t < e.length; ++t) o(t);
		return r
	}
}

function DataTextureLoader(e) {
	Loader.call(this, e)
}

function TextureLoader(e) {
	Loader.call(this, e)
}

function Curve() {
	this.type = "Curve", this.arcLengthDivisions = 200
}
DataTextureLoader.prototype = Object.assign(Object.create(Loader.prototype), {
	constructor: DataTextureLoader,
	load: function(e, t, i, n) {
		const r = this,
			a = new DataTexture,
			s = new FileLoader(this.manager);
		return s.setResponseType("arraybuffer"), s.setRequestHeader(this.requestHeader), s.setPath(this.path), s.setWithCredentials(r.withCredentials), s.load(e, (function(e) {
			const i = r.parse(e);
			i && (void 0 !== i.image ? a.image = i.image : void 0 !== i.data && (a.image.width = i.width, a.image.height = i.height, a.image.data = i.data), a.wrapS = void 0 !== i.wrapS ? i.wrapS : 1001, a.wrapT = void 0 !== i.wrapT ? i.wrapT : 1001, a.magFilter = void 0 !== i.magFilter ? i.magFilter : 1006, a.minFilter = void 0 !== i.minFilter ? i.minFilter : 1006, a.anisotropy = void 0 !== i.anisotropy ? i.anisotropy : 1, void 0 !== i.encoding && (a.encoding = i.encoding), void 0 !== i.flipY && (a.flipY = i.flipY), void 0 !== i.format && (a.format = i.format), void 0 !== i.type && (a.type = i.type), void 0 !== i.mipmaps && (a.mipmaps = i.mipmaps, a.minFilter = 1008), 1 === i.mipmapCount && (a.minFilter = 1006), a.needsUpdate = !0, t && t(a, i))
		}), i, n), a
	}
}), TextureLoader.prototype = Object.assign(Object.create(Loader.prototype), {
	constructor: TextureLoader,
	load: function(e, t, i, n) {
		const r = new Texture$1,
			a = new ImageLoader(this.manager);
		return a.setCrossOrigin(this.crossOrigin), a.setPath(this.path), a.load(e, (function(i) {
			r.image = i;
			const n = e.search(/\.jpe?g($|\?)/i) > 0 || 0 === e.search(/^data\:image\/jpeg/);
			r.format = n ? 1022 : 1023, r.needsUpdate = !0, void 0 !== t && t(r)
		}), i, n), r
	}
}), Object.assign(Curve.prototype, {
	getPoint: function() {
		return console.warn("THREE.Curve: .getPoint() not implemented."), null
	},
	getPointAt: function(e, t) {
		const i = this.getUtoTmapping(e);
		return this.getPoint(i, t)
	},
	getPoints: function(e = 5) {
		const t = [];
		for (let i = 0; i <= e; i++) t.push(this.getPoint(i / e));
		return t
	},
	getSpacedPoints: function(e = 5) {
		const t = [];
		for (let i = 0; i <= e; i++) t.push(this.getPointAt(i / e));
		return t
	},
	getLength: function() {
		const e = this.getLengths();
		return e[e.length - 1]
	},
	getLengths: function(e) {
		if (void 0 === e && (e = this.arcLengthDivisions), this.cacheArcLengths && this.cacheArcLengths.length === e + 1 && !this.needsUpdate) return this.cacheArcLengths;
		this.needsUpdate = !1;
		const t = [];
		let i, n = this.getPoint(0),
			r = 0;
		t.push(0);
		for (let a = 1; a <= e; a++) i = this.getPoint(a / e), r += i.distanceTo(n), t.push(r), n = i;
		return this.cacheArcLengths = t, t
	},
	updateArcLengths: function() {
		this.needsUpdate = !0, this.getLengths()
	},
	getUtoTmapping: function(e, t) {
		const i = this.getLengths();
		let n = 0;
		const r = i.length;
		let a;
		a = t || e * i[r - 1];
		let s, o = 0,
			l = r - 1;
		for (; o <= l;)
			if (n = Math.floor(o + (l - o) / 2), s = i[n] - a, s < 0) o = n + 1;
			else {
				if (!(s > 0)) {
					l = n;
					break
				}
				l = n - 1
			} if (n = l, i[n] === a) return n / (r - 1);
		const c = i[n];
		return (n + (a - c) / (i[n + 1] - c)) / (r - 1)
	},
	getTangent: function(e, t) {
		let i = e - 1e-4,
			n = e + 1e-4;
		i < 0 && (i = 0), n > 1 && (n = 1);
		const r = this.getPoint(i),
			a = this.getPoint(n),
			s = t || (r.isVector2 ? new Vector2 : new Vector3);
		return s.copy(a).sub(r).normalize(), s
	},
	getTangentAt: function(e, t) {
		const i = this.getUtoTmapping(e);
		return this.getTangent(i, t)
	},
	computeFrenetFrames: function(e, t) {
		const i = new Vector3,
			n = [],
			r = [],
			a = [],
			s = new Vector3,
			o = new Matrix4;
		for (let t = 0; t <= e; t++) {
			const i = t / e;
			n[t] = this.getTangentAt(i, new Vector3), n[t].normalize()
		}
		r[0] = new Vector3, a[0] = new Vector3;
		let l = Number.MAX_VALUE;
		const c = Math.abs(n[0].x),
			h = Math.abs(n[0].y),
			u = Math.abs(n[0].z);
		c <= l && (l = c, i.set(1, 0, 0)), h <= l && (l = h, i.set(0, 1, 0)), u <= l && i.set(0, 0, 1), s.crossVectors(n[0], i).normalize(), r[0].crossVectors(n[0], s), a[0].crossVectors(n[0], r[0]);
		for (let t = 1; t <= e; t++) {
			if (r[t] = r[t - 1].clone(), a[t] = a[t - 1].clone(), s.crossVectors(n[t - 1], n[t]), s.length() > Number.EPSILON) {
				s.normalize();
				const e = Math.acos(MathUtils.clamp(n[t - 1].dot(n[t]), -1, 1));
				r[t].applyMatrix4(o.makeRotationAxis(s, e))
			}
			a[t].crossVectors(n[t], r[t])
		}
		if (!0 === t) {
			let t = Math.acos(MathUtils.clamp(r[0].dot(r[e]), -1, 1));
			t /= e, n[0].dot(s.crossVectors(r[0], r[e])) > 0 && (t = -t);
			for (let i = 1; i <= e; i++) r[i].applyMatrix4(o.makeRotationAxis(n[i], t * i)), a[i].crossVectors(n[i], r[i])
		}
		return {
			tangents: n,
			normals: r,
			binormals: a
		}
	},
	clone: function() {
		return (new this.constructor).copy(this)
	},
	copy: function(e) {
		return this.arcLengthDivisions = e.arcLengthDivisions, this
	},
	toJSON: function() {
		const e = {
			metadata: {
				version: 4.5,
				type: "Curve",
				generator: "Curve.toJSON"
			}
		};
		return e.arcLengthDivisions = this.arcLengthDivisions, e.type = this.type, e
	},
	fromJSON: function(e) {
		return this.arcLengthDivisions = e.arcLengthDivisions, this
	}
});
class Light extends Object3D {
	constructor(e, t = 1) {
		super(), this.type = "Light", this.color = new Color(e), this.intensity = t
	}
	copy(e) {
		return super.copy(e), this.color.copy(e.color), this.intensity = e.intensity, this
	}
	toJSON(e) {
		const t = super.toJSON(e);
		return t.object.color = this.color.getHex(), t.object.intensity = this.intensity, void 0 !== this.groundColor && (t.object.groundColor = this.groundColor.getHex()), void 0 !== this.distance && (t.object.distance = this.distance), void 0 !== this.angle && (t.object.angle = this.angle), void 0 !== this.decay && (t.object.decay = this.decay), void 0 !== this.penumbra && (t.object.penumbra = this.penumbra), void 0 !== this.shadow && (t.object.shadow = this.shadow.toJSON()), t
	}
}
Light.prototype.isLight = !0;
const _projScreenMatrix = new Matrix4,
	_lightPositionWorld = new Vector3,
	_lookTarget = new Vector3;
class LightShadow {
	constructor(e) {
		this.camera = e, this.bias = 0, this.normalBias = 0, this.radius = 1, this.mapSize = new Vector2(512, 512), this.map = null, this.mapPass = null, this.matrix = new Matrix4, this.autoUpdate = !0, this.needsUpdate = !1, this._frustum = new Frustum, this._frameExtents = new Vector2(1, 1), this._viewportCount = 1, this._viewports = [new Vector4(0, 0, 1, 1)]
	}
	getViewportCount() {
		return this._viewportCount
	}
	getFrustum() {
		return this._frustum
	}
	updateMatrices(e) {
		const t = this.camera,
			i = this.matrix;
		_lightPositionWorld.setFromMatrixPosition(e.matrixWorld), t.position.copy(_lightPositionWorld), _lookTarget.setFromMatrixPosition(e.target.matrixWorld), t.lookAt(_lookTarget), t.updateMatrixWorld(), _projScreenMatrix.multiplyMatrices(t.projectionMatrix, t.matrixWorldInverse), this._frustum.setFromProjectionMatrix(_projScreenMatrix), i.set(.5, 0, 0, .5, 0, .5, 0, .5, 0, 0, .5, .5, 0, 0, 0, 1), i.multiply(t.projectionMatrix), i.multiply(t.matrixWorldInverse)
	}
	getViewport(e) {
		return this._viewports[e]
	}
	getFrameExtents() {
		return this._frameExtents
	}
	copy(e) {
		return this.camera = e.camera.clone(), this.bias = e.bias, this.radius = e.radius, this.mapSize.copy(e.mapSize), this
	}
	clone() {
		return (new this.constructor).copy(this)
	}
	toJSON() {
		const e = {};
		return 0 !== this.bias && (e.bias = this.bias), 0 !== this.normalBias && (e.normalBias = this.normalBias), 1 !== this.radius && (e.radius = this.radius), 512 === this.mapSize.x && 512 === this.mapSize.y || (e.mapSize = this.mapSize.toArray()), e.camera = this.camera.toJSON(!1).object, delete e.camera.matrix, e
	}
}
class SpotLightShadow extends LightShadow {
	constructor() {
		super(new PerspectiveCamera(50, 1, .5, 500)), this.focus = 1
	}
	updateMatrices(e) {
		const t = this.camera,
			i = 2 * MathUtils.RAD2DEG * e.angle * this.focus,
			n = this.mapSize.width / this.mapSize.height,
			r = e.distance || t.far;
		i === t.fov && n === t.aspect && r === t.far || (t.fov = i, t.aspect = n, t.far = r, t.updateProjectionMatrix()), super.updateMatrices(e)
	}
}
SpotLightShadow.prototype.isSpotLightShadow = !0;
class SpotLight extends Light {
	constructor(e, t, i = 0, n = Math.PI / 3, r = 0, a = 1) {
		super(e, t), this.type = "SpotLight", this.position.copy(Object3D.DefaultUp), this.updateMatrix(), this.target = new Object3D, this.distance = i, this.angle = n, this.penumbra = r, this.decay = a, this.shadow = new SpotLightShadow
	}
	get power() {
		return this.intensity * Math.PI
	}
	set power(e) {
		this.intensity = e / Math.PI
	}
	copy(e) {
		return super.copy(e), this.distance = e.distance, this.angle = e.angle, this.penumbra = e.penumbra, this.decay = e.decay, this.target = e.target.clone(), this.shadow = e.shadow.clone(), this
	}
}
SpotLight.prototype.isSpotLight = !0;
const _projScreenMatrix$1 = new Matrix4,
	_lightPositionWorld$1 = new Vector3,
	_lookTarget$1 = new Vector3;
class PointLightShadow extends LightShadow {
	constructor() {
		super(new PerspectiveCamera(90, 1, .5, 500)), this._frameExtents = new Vector2(4, 2), this._viewportCount = 6, this._viewports = [new Vector4(2, 1, 1, 1), new Vector4(0, 1, 1, 1), new Vector4(3, 1, 1, 1), new Vector4(1, 1, 1, 1), new Vector4(3, 0, 1, 1), new Vector4(1, 0, 1, 1)], this._cubeDirections = [new Vector3(1, 0, 0), new Vector3(-1, 0, 0), new Vector3(0, 0, 1), new Vector3(0, 0, -1), new Vector3(0, 1, 0), new Vector3(0, -1, 0)], this._cubeUps = [new Vector3(0, 1, 0), new Vector3(0, 1, 0), new Vector3(0, 1, 0), new Vector3(0, 1, 0), new Vector3(0, 0, 1), new Vector3(0, 0, -1)]
	}
	updateMatrices(e, t = 0) {
		const i = this.camera,
			n = this.matrix;
		_lightPositionWorld$1.setFromMatrixPosition(e.matrixWorld), i.position.copy(_lightPositionWorld$1), _lookTarget$1.copy(i.position), _lookTarget$1.add(this._cubeDirections[t]), i.up.copy(this._cubeUps[t]), i.lookAt(_lookTarget$1), i.updateMatrixWorld(), n.makeTranslation(-_lightPositionWorld$1.x, -_lightPositionWorld$1.y, -_lightPositionWorld$1.z), _projScreenMatrix$1.multiplyMatrices(i.projectionMatrix, i.matrixWorldInverse), this._frustum.setFromProjectionMatrix(_projScreenMatrix$1)
	}
}
PointLightShadow.prototype.isPointLightShadow = !0;
class PointLight extends Light {
	constructor(e, t, i = 0, n = 1) {
		super(e, t), this.type = "PointLight", this.distance = i, this.decay = n, this.shadow = new PointLightShadow
	}
	get power() {
		return 4 * this.intensity * Math.PI
	}
	set power(e) {
		this.intensity = e / (4 * Math.PI)
	}
	copy(e) {
		return super.copy(e), this.distance = e.distance, this.decay = e.decay, this.shadow = e.shadow.clone(), this
	}
}
PointLight.prototype.isPointLight = !0;
class OrthographicCamera extends Camera {
	constructor(e = -1, t = 1, i = 1, n = -1, r = .1, a = 2e3) {
		super(), this.type = "OrthographicCamera", this.zoom = 1, this.view = null, this.left = e, this.right = t, this.top = i, this.bottom = n, this.near = r, this.far = a, this.updateProjectionMatrix()
	}
	copy(e, t) {
		return super.copy(e, t), this.left = e.left, this.right = e.right, this.top = e.top, this.bottom = e.bottom, this.near = e.near, this.far = e.far, this.zoom = e.zoom, this.view = null === e.view ? null : Object.assign({}, e.view), this
	}
	setViewOffset(e, t, i, n, r, a) {
		null === this.view && (this.view = {
			enabled: !0,
			fullWidth: 1,
			fullHeight: 1,
			offsetX: 0,
			offsetY: 0,
			width: 1,
			height: 1
		}), this.view.enabled = !0, this.view.fullWidth = e, this.view.fullHeight = t, this.view.offsetX = i, this.view.offsetY = n, this.view.width = r, this.view.height = a, this.updateProjectionMatrix()
	}
	clearViewOffset() {
		null !== this.view && (this.view.enabled = !1), this.updateProjectionMatrix()
	}
	updateProjectionMatrix() {
		const e = (this.right - this.left) / (2 * this.zoom),
			t = (this.top - this.bottom) / (2 * this.zoom),
			i = (this.right + this.left) / 2,
			n = (this.top + this.bottom) / 2;
		let r = i - e,
			a = i + e,
			s = n + t,
			o = n - t;
		if (null !== this.view && this.view.enabled) {
			const e = (this.right - this.left) / this.view.fullWidth / this.zoom,
				t = (this.top - this.bottom) / this.view.fullHeight / this.zoom;
			r += e * this.view.offsetX, a = r + e * this.view.width, s -= t * this.view.offsetY, o = s - t * this.view.height
		}
		this.projectionMatrix.makeOrthographic(r, a, s, o, this.near, this.far), this.projectionMatrixInverse.copy(this.projectionMatrix).invert()
	}
	toJSON(e) {
		const t = Object3D.prototype.toJSON.call(this, e);
		return t.object.zoom = this.zoom, t.object.left = this.left, t.object.right = this.right, t.object.top = this.top, t.object.bottom = this.bottom, t.object.near = this.near, t.object.far = this.far, null !== this.view && (t.object.view = Object.assign({}, this.view)), t
	}
}
OrthographicCamera.prototype.isOrthographicCamera = !0;
class DirectionalLightShadow extends LightShadow {
	constructor() {
		super(new OrthographicCamera(-5, 5, 5, -5, .5, 500))
	}
}
DirectionalLightShadow.prototype.isDirectionalLightShadow = !0;
class DirectionalLight extends Light {
	constructor(e, t) {
		super(e, t), this.type = "DirectionalLight", this.position.copy(Object3D.DefaultUp), this.updateMatrix(), this.target = new Object3D, this.shadow = new DirectionalLightShadow
	}
	copy(e) {
		return super.copy(e), this.target = e.target.clone(), this.shadow = e.shadow.clone(), this
	}
}
DirectionalLight.prototype.isDirectionalLight = !0;
const LoaderUtils = {
	decodeText: function(e) {
		if ("undefined" != typeof TextDecoder) return (new TextDecoder).decode(e);
		let t = "";
		for (let i = 0, n = e.length; i < n; i++) t += String.fromCharCode(e[i]);
		try {
			return decodeURIComponent(escape(t))
		} catch (e) {
			return t
		}
	},
	extractUrlBase: function(e) {
		const t = e.lastIndexOf("/");
		return -1 === t ? "./" : e.substr(0, t + 1)
	}
};

function InstancedBufferGeometry() {
	BufferGeometry.call(this), this.type = "InstancedBufferGeometry", this.instanceCount = 1 / 0
}

function InstancedBufferAttribute(e, t, i, n) {
	"number" == typeof i && (n = i, i = !1, console.error("THREE.InstancedBufferAttribute: The constructor now expects normalized as the third argument.")), BufferAttribute.call(this, e, t, i), this.meshPerAttribute = n || 1
}

function ImageBitmapLoader(e) {
	"undefined" == typeof createImageBitmap && console.warn("THREE.ImageBitmapLoader: createImageBitmap() not supported."), "undefined" == typeof fetch && console.warn("THREE.ImageBitmapLoader: fetch() not supported."), Loader.call(this, e), this.options = {
		premultiplyAlpha: "none"
	}
}
let _context;
InstancedBufferGeometry.prototype = Object.assign(Object.create(BufferGeometry.prototype), {
	constructor: InstancedBufferGeometry,
	isInstancedBufferGeometry: !0,
	copy: function(e) {
		return BufferGeometry.prototype.copy.call(this, e), this.instanceCount = e.instanceCount, this
	},
	clone: function() {
		return (new this.constructor).copy(this)
	},
	toJSON: function() {
		const e = BufferGeometry.prototype.toJSON.call(this);
		return e.instanceCount = this.instanceCount, e.isInstancedBufferGeometry = !0, e
	}
}), InstancedBufferAttribute.prototype = Object.assign(Object.create(BufferAttribute.prototype), {
	constructor: InstancedBufferAttribute,
	isInstancedBufferAttribute: !0,
	copy: function(e) {
		return BufferAttribute.prototype.copy.call(this, e), this.meshPerAttribute = e.meshPerAttribute, this
	},
	toJSON: function() {
		const e = BufferAttribute.prototype.toJSON.call(this);
		return e.meshPerAttribute = this.meshPerAttribute, e.isInstancedBufferAttribute = !0, e
	}
}), ImageBitmapLoader.prototype = Object.assign(Object.create(Loader.prototype), {
	constructor: ImageBitmapLoader,
	isImageBitmapLoader: !0,
	setOptions: function(e) {
		return this.options = e, this
	},
	load: function(e, t, i, n) {
		void 0 === e && (e = ""), void 0 !== this.path && (e = this.path + e), e = this.manager.resolveURL(e);
		const r = this,
			a = Cache.get(e);
		if (void 0 !== a) return r.manager.itemStart(e), setTimeout((function() {
			t && t(a), r.manager.itemEnd(e)
		}), 0), a;
		const s = {};
		s.credentials = "anonymous" === this.crossOrigin ? "same-origin" : "include", s.headers = this.requestHeader, fetch(e, s).then((function(e) {
			return e.blob()
		})).then((function(e) {
			return createImageBitmap(e, Object.assign(r.options, {
				colorSpaceConversion: "none"
			}))
		})).then((function(i) {
			Cache.add(e, i), t && t(i), r.manager.itemEnd(e)
		})).catch((function(t) {
			n && n(t), r.manager.itemError(e), r.manager.itemEnd(e)
		})), r.manager.itemStart(e)
	}
});
const AudioContext = {
	getContext: function() {
		return void 0 === _context && (_context = new(window.AudioContext || window.webkitAudioContext)), _context
	},
	setContext: function(e) {
		_context = e
	}
};
class AudioLoader extends Loader {
	constructor(e) {
		super(e)
	}
	load(e, t, i, n) {
		const r = this,
			a = new FileLoader(this.manager);
		a.setResponseType("arraybuffer"), a.setPath(this.path), a.setRequestHeader(this.requestHeader), a.setWithCredentials(this.withCredentials), a.load(e, (function(i) {
			try {
				const e = i.slice(0);
				AudioContext.getContext().decodeAudioData(e, (function(e) {
					t(e)
				}))
			} catch (t) {
				n ? n(t) : console.error(t), r.manager.itemError(e)
			}
		}), i, n)
	}
}
new Matrix4, new Matrix4;
class Audio extends Object3D {
	constructor(e) {
		super(), this.type = "Audio", this.listener = e, this.context = e.context, this.gain = this.context.createGain(), this.gain.connect(e.getInput()), this.autoplay = !1, this.buffer = null, this.detune = 0, this.loop = !1, this.loopStart = 0, this.loopEnd = 0, this.offset = 0, this.duration = void 0, this.playbackRate = 1, this.isPlaying = !1, this.hasPlaybackControl = !0, this.source = null, this.sourceType = "empty", this._startedAt = 0, this._progress = 0, this._connected = !1, this.filters = []
	}
	getOutput() {
		return this.gain
	}
	setNodeSource(e) {
		return this.hasPlaybackControl = !1, this.sourceType = "audioNode", this.source = e, this.connect(), this
	}
	setMediaElementSource(e) {
		return this.hasPlaybackControl = !1, this.sourceType = "mediaNode", this.source = this.context.createMediaElementSource(e), this.connect(), this
	}
	setMediaStreamSource(e) {
		return this.hasPlaybackControl = !1, this.sourceType = "mediaStreamNode", this.source = this.context.createMediaStreamSource(e), this.connect(), this
	}
	setBuffer(e) {
		return this.buffer = e, this.sourceType = "buffer", this.autoplay && this.play(), this
	}
	play(e = 0) {
		if (!0 === this.isPlaying) return void console.warn("THREE.Audio: Audio is already playing.");
		if (!1 === this.hasPlaybackControl) return void console.warn("THREE.Audio: this Audio has no playback control.");
		this._startedAt = this.context.currentTime + e;
		const t = this.context.createBufferSource();
		return t.buffer = this.buffer, t.loop = this.loop, t.loopStart = this.loopStart, t.loopEnd = this.loopEnd, t.onended = this.onEnded.bind(this), t.start(this._startedAt, this._progress + this.offset, this.duration), this.isPlaying = !0, this.source = t, this.setDetune(this.detune), this.setPlaybackRate(this.playbackRate), this.connect()
	}
	pause() {
		if (!1 !== this.hasPlaybackControl) return !0 === this.isPlaying && (this._progress += Math.max(this.context.currentTime - this._startedAt, 0) * this.playbackRate, !0 === this.loop && (this._progress = this._progress % (this.duration || this.buffer.duration)), this.source.stop(), this.source.onended = null, this.isPlaying = !1), this;
		console.warn("THREE.Audio: this Audio has no playback control.")
	}
	stop() {
		if (!1 !== this.hasPlaybackControl) return this._progress = 0, this.source.stop(), this.source.onended = null, this.isPlaying = !1, this;
		console.warn("THREE.Audio: this Audio has no playback control.")
	}
	connect() {
		if (this.filters.length > 0) {
			this.source.connect(this.filters[0]);
			for (let e = 1, t = this.filters.length; e < t; e++) this.filters[e - 1].connect(this.filters[e]);
			this.filters[this.filters.length - 1].connect(this.getOutput())
		} else this.source.connect(this.getOutput());
		return this._connected = !0, this
	}
	disconnect() {
		if (this.filters.length > 0) {
			this.source.disconnect(this.filters[0]);
			for (let e = 1, t = this.filters.length; e < t; e++) this.filters[e - 1].disconnect(this.filters[e]);
			this.filters[this.filters.length - 1].disconnect(this.getOutput())
		} else this.source.disconnect(this.getOutput());
		return this._connected = !1, this
	}
	getFilters() {
		return this.filters
	}
	setFilters(e) {
		return e || (e = []), !0 === this._connected ? (this.disconnect(), this.filters = e.slice(), this.connect()) : this.filters = e.slice(), this
	}
	setDetune(e) {
		if (this.detune = e, void 0 !== this.source.detune) return !0 === this.isPlaying && this.source.detune.setTargetAtTime(this.detune, this.context.currentTime, .01), this
	}
	getDetune() {
		return this.detune
	}
	getFilter() {
		return this.getFilters()[0]
	}
	setFilter(e) {
		return this.setFilters(e ? [e] : [])
	}
	setPlaybackRate(e) {
		if (!1 !== this.hasPlaybackControl) return this.playbackRate = e, !0 === this.isPlaying && this.source.playbackRate.setTargetAtTime(this.playbackRate, this.context.currentTime, .01), this;
		console.warn("THREE.Audio: this Audio has no playback control.")
	}
	getPlaybackRate() {
		return this.playbackRate
	}
	onEnded() {
		this.isPlaying = !1
	}
	getLoop() {
		return !1 === this.hasPlaybackControl ? (console.warn("THREE.Audio: this Audio has no playback control."), !1) : this.loop
	}
	setLoop(e) {
		if (!1 !== this.hasPlaybackControl) return this.loop = e, !0 === this.isPlaying && (this.source.loop = this.loop), this;
		console.warn("THREE.Audio: this Audio has no playback control.")
	}
	setLoopStart(e) {
		return this.loopStart = e, this
	}
	setLoopEnd(e) {
		return this.loopEnd = e, this
	}
	getVolume() {
		return this.gain.gain.value
	}
	setVolume(e) {
		return this.gain.gain.setTargetAtTime(e, this.context.currentTime, .01), this
	}
}
class PropertyMixer {
	constructor(e, t, i) {
		let n, r, a;
		switch (this.binding = e, this.valueSize = i, t) {
			case "quaternion":
				n = this._slerp, r = this._slerpAdditive, a = this._setAdditiveIdentityQuaternion, this.buffer = new Float64Array(6 * i), this._workIndex = 5;
				break;
			case "string":
			case "bool":
				n = this._select, r = this._select, a = this._setAdditiveIdentityOther, this.buffer = new Array(5 * i);
				break;
			default:
				n = this._lerp, r = this._lerpAdditive, a = this._setAdditiveIdentityNumeric, this.buffer = new Float64Array(5 * i)
		}
		this._mixBufferRegion = n, this._mixBufferRegionAdditive = r, this._setIdentity = a, this._origIndex = 3, this._addIndex = 4, this.cumulativeWeight = 0, this.cumulativeWeightAdditive = 0, this.useCount = 0, this.referenceCount = 0
	}
	accumulate(e, t) {
		const i = this.buffer,
			n = this.valueSize,
			r = e * n + n;
		let a = this.cumulativeWeight;
		if (0 === a) {
			for (let e = 0; e !== n; ++e) i[r + e] = i[e];
			a = t
		} else {
			a += t;
			const e = t / a;
			this._mixBufferRegion(i, r, 0, e, n)
		}
		this.cumulativeWeight = a
	}
	accumulateAdditive(e) {
		const t = this.buffer,
			i = this.valueSize,
			n = i * this._addIndex;
		0 === this.cumulativeWeightAdditive && this._setIdentity(), this._mixBufferRegionAdditive(t, n, 0, e, i), this.cumulativeWeightAdditive += e
	}
	apply(e) {
		const t = this.valueSize,
			i = this.buffer,
			n = e * t + t,
			r = this.cumulativeWeight,
			a = this.cumulativeWeightAdditive,
			s = this.binding;
		if (this.cumulativeWeight = 0, this.cumulativeWeightAdditive = 0, r < 1) {
			const e = t * this._origIndex;
			this._mixBufferRegion(i, n, e, 1 - r, t)
		}
		a > 0 && this._mixBufferRegionAdditive(i, n, this._addIndex * t, 1, t);
		for (let e = t, r = t + t; e !== r; ++e)
			if (i[e] !== i[e + t]) {
				s.setValue(i, n);
				break
			}
	}
	saveOriginalState() {
		const e = this.binding,
			t = this.buffer,
			i = this.valueSize,
			n = i * this._origIndex;
		e.getValue(t, n);
		for (let e = i, r = n; e !== r; ++e) t[e] = t[n + e % i];
		this._setIdentity(), this.cumulativeWeight = 0, this.cumulativeWeightAdditive = 0
	}
	restoreOriginalState() {
		const e = 3 * this.valueSize;
		this.binding.setValue(this.buffer, e)
	}
	_setAdditiveIdentityNumeric() {
		const e = this._addIndex * this.valueSize,
			t = e + this.valueSize;
		for (let i = e; i < t; i++) this.buffer[i] = 0
	}
	_setAdditiveIdentityQuaternion() {
		this._setAdditiveIdentityNumeric(), this.buffer[this._addIndex * this.valueSize + 3] = 1
	}
	_setAdditiveIdentityOther() {
		const e = this._origIndex * this.valueSize,
			t = this._addIndex * this.valueSize;
		for (let i = 0; i < this.valueSize; i++) this.buffer[t + i] = this.buffer[e + i]
	}
	_select(e, t, i, n, r) {
		if (n >= .5)
			for (let n = 0; n !== r; ++n) e[t + n] = e[i + n]
	}
	_slerp(e, t, i, n) {
		Quaternion.slerpFlat(e, t, e, t, e, i, n)
	}
	_slerpAdditive(e, t, i, n, r) {
		const a = this._workIndex * r;
		Quaternion.multiplyQuaternionsFlat(e, a, e, t, e, i), Quaternion.slerpFlat(e, t, e, t, e, a, n)
	}
	_lerp(e, t, i, n, r) {
		const a = 1 - n;
		for (let s = 0; s !== r; ++s) {
			const r = t + s;
			e[r] = e[r] * a + e[i + s] * n
		}
	}
	_lerpAdditive(e, t, i, n, r) {
		for (let a = 0; a !== r; ++a) {
			const r = t + a;
			e[r] = e[r] + e[i + a] * n
		}
	}
}
const _RESERVED_CHARS_RE = "\\[\\]\\.:\\/",
	_reservedRe = new RegExp("[\\[\\]\\.:\\/]", "g"),
	_wordChar = "[^\\[\\]\\.:\\/]",
	_wordCharOrDot = "[^" + "\\[\\]\\.:\\/".replace("\\.", "") + "]",
	_directoryRe = /((?:WC+[\/:])*)/.source.replace("WC", _wordChar),
	_nodeRe = /(WCOD+)?/.source.replace("WCOD", _wordCharOrDot),
	_objectRe = /(?:\.(WC+)(?:\[(.+)\])?)?/.source.replace("WC", _wordChar),
	_propertyRe = /\.(WC+)(?:\[(.+)\])?/.source.replace("WC", _wordChar),
	_trackRe = new RegExp("^" + _directoryRe + _nodeRe + _objectRe + _propertyRe + "$"),
	_supportedObjectNames = ["material", "materials", "bones"];

function Composite(e, t, i) {
	const n = i || PropertyBinding.parseTrackName(t);
	this._targetGroup = e, this._bindings = e.subscribe_(t, n)
}

function PropertyBinding(e, t, i) {
	this.path = t, this.parsedPath = i || PropertyBinding.parseTrackName(t), this.node = PropertyBinding.findNode(e, this.parsedPath.nodeName) || e, this.rootNode = e
}
Object.assign(Composite.prototype, {
	getValue: function(e, t) {
		this.bind();
		const i = this._targetGroup.nCachedObjects_,
			n = this._bindings[i];
		void 0 !== n && n.getValue(e, t)
	},
	setValue: function(e, t) {
		const i = this._bindings;
		for (let n = this._targetGroup.nCachedObjects_, r = i.length; n !== r; ++n) i[n].setValue(e, t)
	},
	bind: function() {
		const e = this._bindings;
		for (let t = this._targetGroup.nCachedObjects_, i = e.length; t !== i; ++t) e[t].bind()
	},
	unbind: function() {
		const e = this._bindings;
		for (let t = this._targetGroup.nCachedObjects_, i = e.length; t !== i; ++t) e[t].unbind()
	}
}), Object.assign(PropertyBinding, {
	Composite: Composite,
	create: function(e, t, i) {
		return e && e.isAnimationObjectGroup ? new PropertyBinding.Composite(e, t, i) : new PropertyBinding(e, t, i)
	},
	sanitizeNodeName: function(e) {
		return e.replace(/\s/g, "_").replace(_reservedRe, "")
	},
	parseTrackName: function(e) {
		const t = _trackRe.exec(e);
		if (!t) throw new Error("PropertyBinding: Cannot parse trackName: " + e);
		const i = {
				nodeName: t[2],
				objectName: t[3],
				objectIndex: t[4],
				propertyName: t[5],
				propertyIndex: t[6]
			},
			n = i.nodeName && i.nodeName.lastIndexOf(".");
		if (void 0 !== n && -1 !== n) {
			const e = i.nodeName.substring(n + 1); - 1 !== _supportedObjectNames.indexOf(e) && (i.nodeName = i.nodeName.substring(0, n), i.objectName = e)
		}
		if (null === i.propertyName || 0 === i.propertyName.length) throw new Error("PropertyBinding: can not parse propertyName from trackName: " + e);
		return i
	},
	findNode: function(e, t) {
		if (!t || "" === t || "." === t || -1 === t || t === e.name || t === e.uuid) return e;
		if (e.skeleton) {
			const i = e.skeleton.getBoneByName(t);
			if (void 0 !== i) return i
		}
		if (e.children) {
			const i = function(e) {
					for (let n = 0; n < e.length; n++) {
						const r = e[n];
						if (r.name === t || r.uuid === t) return r;
						const a = i(r.children);
						if (a) return a
					}
					return null
				},
				n = i(e.children);
			if (n) return n
		}
		return null
	}
}), Object.assign(PropertyBinding.prototype, {
	_getValue_unavailable: function() {},
	_setValue_unavailable: function() {},
	BindingType: {
		Direct: 0,
		EntireArray: 1,
		ArrayElement: 2,
		HasFromToArray: 3
	},
	Versioning: {
		None: 0,
		NeedsUpdate: 1,
		MatrixWorldNeedsUpdate: 2
	},
	GetterByBindingType: [function(e, t) {
		e[t] = this.node[this.propertyName]
	}, function(e, t) {
		const i = this.resolvedProperty;
		for (let n = 0, r = i.length; n !== r; ++n) e[t++] = i[n]
	}, function(e, t) {
		e[t] = this.resolvedProperty[this.propertyIndex]
	}, function(e, t) {
		this.resolvedProperty.toArray(e, t)
	}],
	SetterByBindingTypeAndVersioning: [
		[function(e, t) {
			this.targetObject[this.propertyName] = e[t]
		}, function(e, t) {
			this.targetObject[this.propertyName] = e[t], this.targetObject.needsUpdate = !0
		}, function(e, t) {
			this.targetObject[this.propertyName] = e[t], this.targetObject.matrixWorldNeedsUpdate = !0
		}],
		[function(e, t) {
			const i = this.resolvedProperty;
			for (let n = 0, r = i.length; n !== r; ++n) i[n] = e[t++]
		}, function(e, t) {
			const i = this.resolvedProperty;
			for (let n = 0, r = i.length; n !== r; ++n) i[n] = e[t++];
			this.targetObject.needsUpdate = !0
		}, function(e, t) {
			const i = this.resolvedProperty;
			for (let n = 0, r = i.length; n !== r; ++n) i[n] = e[t++];
			this.targetObject.matrixWorldNeedsUpdate = !0
		}],
		[function(e, t) {
			this.resolvedProperty[this.propertyIndex] = e[t]
		}, function(e, t) {
			this.resolvedProperty[this.propertyIndex] = e[t], this.targetObject.needsUpdate = !0
		}, function(e, t) {
			this.resolvedProperty[this.propertyIndex] = e[t], this.targetObject.matrixWorldNeedsUpdate = !0
		}],
		[function(e, t) {
			this.resolvedProperty.fromArray(e, t)
		}, function(e, t) {
			this.resolvedProperty.fromArray(e, t), this.targetObject.needsUpdate = !0
		}, function(e, t) {
			this.resolvedProperty.fromArray(e, t), this.targetObject.matrixWorldNeedsUpdate = !0
		}]
	],
	getValue: function(e, t) {
		this.bind(), this.getValue(e, t)
	},
	setValue: function(e, t) {
		this.bind(), this.setValue(e, t)
	},
	bind: function() {
		let e = this.node;
		const t = this.parsedPath,
			i = t.objectName,
			n = t.propertyName;
		let r = t.propertyIndex;
		if (e || (e = PropertyBinding.findNode(this.rootNode, t.nodeName) || this.rootNode, this.node = e), this.getValue = this._getValue_unavailable, this.setValue = this._setValue_unavailable, !e) return void console.error("THREE.PropertyBinding: Trying to update node for track: " + this.path + " but it wasn't found.");
		if (i) {
			let n = t.objectIndex;
			switch (i) {
				case "materials":
					if (!e.material) return void console.error("THREE.PropertyBinding: Can not bind to material as node does not have a material.", this);
					if (!e.material.materials) return void console.error("THREE.PropertyBinding: Can not bind to material.materials as node.material does not have a materials array.", this);
					e = e.material.materials;
					break;
				case "bones":
					if (!e.skeleton) return void console.error("THREE.PropertyBinding: Can not bind to bones as node does not have a skeleton.", this);
					e = e.skeleton.bones;
					for (let t = 0; t < e.length; t++)
						if (e[t].name === n) {
							n = t;
							break
						} break;
				default:
					if (void 0 === e[i]) return void console.error("THREE.PropertyBinding: Can not bind to objectName of node undefined.", this);
					e = e[i]
			}
			if (void 0 !== n) {
				if (void 0 === e[n]) return void console.error("THREE.PropertyBinding: Trying to bind to objectIndex of objectName, but is undefined.", this, e);
				e = e[n]
			}
		}
		const a = e[n];
		if (void 0 === a) {
			const i = t.nodeName;
			return void console.error("THREE.PropertyBinding: Trying to update property for track: " + i + "." + n + " but it wasn't found.", e)
		}
		let s = this.Versioning.None;
		this.targetObject = e, void 0 !== e.needsUpdate ? s = this.Versioning.NeedsUpdate : void 0 !== e.matrixWorldNeedsUpdate && (s = this.Versioning.MatrixWorldNeedsUpdate);
		let o = this.BindingType.Direct;
		if (void 0 !== r) {
			if ("morphTargetInfluences" === n) {
				if (!e.geometry) return void console.error("THREE.PropertyBinding: Can not bind to morphTargetInfluences because node does not have a geometry.", this);
				if (!e.geometry.isBufferGeometry) return void console.error("THREE.PropertyBinding: Can not bind to morphTargetInfluences on THREE.Geometry. Use THREE.BufferGeometry instead.", this);
				if (!e.geometry.morphAttributes) return void console.error("THREE.PropertyBinding: Can not bind to morphTargetInfluences because node does not have a geometry.morphAttributes.", this);
				void 0 !== e.morphTargetDictionary[r] && (r = e.morphTargetDictionary[r])
			}
			o = this.BindingType.ArrayElement, this.resolvedProperty = a, this.propertyIndex = r
		} else void 0 !== a.fromArray && void 0 !== a.toArray ? (o = this.BindingType.HasFromToArray, this.resolvedProperty = a) : Array.isArray(a) ? (o = this.BindingType.EntireArray, this.resolvedProperty = a) : this.propertyName = n;
		this.getValue = this.GetterByBindingType[o], this.setValue = this.SetterByBindingTypeAndVersioning[o][s]
	},
	unbind: function() {
		this.node = null, this.getValue = this._getValue_unbound, this.setValue = this._setValue_unbound
	}
}), Object.assign(PropertyBinding.prototype, {
	_getValue_unbound: PropertyBinding.prototype.getValue,
	_setValue_unbound: PropertyBinding.prototype.setValue
});
class AnimationAction {
	constructor(e, t, i = null, n = t.blendMode) {
		this._mixer = e, this._clip = t, this._localRoot = i, this.blendMode = n;
		const r = t.tracks,
			a = r.length,
			s = new Array(a),
			o = {
				endingStart: 2400,
				endingEnd: 2400
			};
		for (let e = 0; e !== a; ++e) {
			const t = r[e].createInterpolant(null);
			s[e] = t, t.settings = o
		}
		this._interpolantSettings = o, this._interpolants = s, this._propertyBindings = new Array(a), this._cacheIndex = null, this._byClipCacheIndex = null, this._timeScaleInterpolant = null, this._weightInterpolant = null, this.loop = 2201, this._loopCount = -1, this._startTime = null, this.time = 0, this.timeScale = 1, this._effectiveTimeScale = 1, this.weight = 1, this._effectiveWeight = 1, this.repetitions = 1 / 0, this.paused = !1, this.enabled = !0, this.clampWhenFinished = !1, this.zeroSlopeAtStart = !0, this.zeroSlopeAtEnd = !0
	}
	play() {
		return this._mixer._activateAction(this), this
	}
	stop() {
		return this._mixer._deactivateAction(this), this.reset()
	}
	reset() {
		return this.paused = !1, this.enabled = !0, this.time = 0, this._loopCount = -1, this._startTime = null, this.stopFading().stopWarping()
	}
	isRunning() {
		return this.enabled && !this.paused && 0 !== this.timeScale && null === this._startTime && this._mixer._isActiveAction(this)
	}
	isScheduled() {
		return this._mixer._isActiveAction(this)
	}
	startAt(e) {
		return this._startTime = e, this
	}
	setLoop(e, t) {
		return this.loop = e, this.repetitions = t, this
	}
	setEffectiveWeight(e) {
		return this.weight = e, this._effectiveWeight = this.enabled ? e : 0, this.stopFading()
	}
	getEffectiveWeight() {
		return this._effectiveWeight
	}
	fadeIn(e) {
		return this._scheduleFading(e, 0, 1)
	}
	fadeOut(e) {
		return this._scheduleFading(e, 1, 0)
	}
	crossFadeFrom(e, t, i) {
		if (e.fadeOut(t), this.fadeIn(t), i) {
			const i = this._clip.duration,
				n = e._clip.duration,
				r = n / i,
				a = i / n;
			e.warp(1, r, t), this.warp(a, 1, t)
		}
		return this
	}
	crossFadeTo(e, t, i) {
		return e.crossFadeFrom(this, t, i)
	}
	stopFading() {
		const e = this._weightInterpolant;
		return null !== e && (this._weightInterpolant = null, this._mixer._takeBackControlInterpolant(e)), this
	}
	setEffectiveTimeScale(e) {
		return this.timeScale = e, this._effectiveTimeScale = this.paused ? 0 : e, this.stopWarping()
	}
	getEffectiveTimeScale() {
		return this._effectiveTimeScale
	}
	setDuration(e) {
		return this.timeScale = this._clip.duration / e, this.stopWarping()
	}
	syncWith(e) {
		return this.time = e.time, this.timeScale = e.timeScale, this.stopWarping()
	}
	halt(e) {
		return this.warp(this._effectiveTimeScale, 0, e)
	}
	warp(e, t, i) {
		const n = this._mixer,
			r = n.time,
			a = this.timeScale;
		let s = this._timeScaleInterpolant;
		null === s && (s = n._lendControlInterpolant(), this._timeScaleInterpolant = s);
		const o = s.parameterPositions,
			l = s.sampleValues;
		return o[0] = r, o[1] = r + i, l[0] = e / a, l[1] = t / a, this
	}
	stopWarping() {
		const e = this._timeScaleInterpolant;
		return null !== e && (this._timeScaleInterpolant = null, this._mixer._takeBackControlInterpolant(e)), this
	}
	getMixer() {
		return this._mixer
	}
	getClip() {
		return this._clip
	}
	getRoot() {
		return this._localRoot || this._mixer._root
	}
	_update(e, t, i, n) {
		if (!this.enabled) return void this._updateWeight(e);
		const r = this._startTime;
		if (null !== r) {
			const n = (e - r) * i;
			if (n < 0 || 0 === i) return;
			this._startTime = null, t = i * n
		}
		t *= this._updateTimeScale(e);
		const a = this._updateTime(t),
			s = this._updateWeight(e);
		if (s > 0) {
			const e = this._interpolants,
				t = this._propertyBindings;
			switch (this.blendMode) {
				case 2501:
					for (let i = 0, n = e.length; i !== n; ++i) e[i].evaluate(a), t[i].accumulateAdditive(s);
					break;
				case 2500:
				default:
					for (let i = 0, r = e.length; i !== r; ++i) e[i].evaluate(a), t[i].accumulate(n, s)
			}
		}
	}
	_updateWeight(e) {
		let t = 0;
		if (this.enabled) {
			t = this.weight;
			const i = this._weightInterpolant;
			if (null !== i) {
				const n = i.evaluate(e)[0];
				t *= n, e > i.parameterPositions[1] && (this.stopFading(), 0 === n && (this.enabled = !1))
			}
		}
		return this._effectiveWeight = t, t
	}
	_updateTimeScale(e) {
		let t = 0;
		if (!this.paused) {
			t = this.timeScale;
			const i = this._timeScaleInterpolant;
			if (null !== i) {
				t *= i.evaluate(e)[0], e > i.parameterPositions[1] && (this.stopWarping(), 0 === t ? this.paused = !0 : this.timeScale = t)
			}
		}
		return this._effectiveTimeScale = t, t
	}
	_updateTime(e) {
		const t = this._clip.duration,
			i = this.loop;
		let n = this.time + e,
			r = this._loopCount;
		const a = 2202 === i;
		if (0 === e) return -1 === r ? n : a && 1 == (1 & r) ? t - n : n;
		if (2200 === i) {
			-1 === r && (this._loopCount = 0, this._setEndings(!0, !0, !1));
			e: {
				if (n >= t) n = t;
				else {
					if (!(n < 0)) {
						this.time = n;
						break e
					}
					n = 0
				}
				this.clampWhenFinished ? this.paused = !0 : this.enabled = !1,
				this.time = n,
				this._mixer.dispatchEvent({
					type: "finished",
					action: this,
					direction: e < 0 ? -1 : 1
				})
			}
		} else {
			if (-1 === r && (e >= 0 ? (r = 0, this._setEndings(!0, 0 === this.repetitions, a)) : this._setEndings(0 === this.repetitions, !0, a)), n >= t || n < 0) {
				const i = Math.floor(n / t);
				n -= t * i, r += Math.abs(i);
				const s = this.repetitions - r;
				if (s <= 0) this.clampWhenFinished ? this.paused = !0 : this.enabled = !1, n = e > 0 ? t : 0, this.time = n, this._mixer.dispatchEvent({
					type: "finished",
					action: this,
					direction: e > 0 ? 1 : -1
				});
				else {
					if (1 === s) {
						const t = e < 0;
						this._setEndings(t, !t, a)
					} else this._setEndings(!1, !1, a);
					this._loopCount = r, this.time = n, this._mixer.dispatchEvent({
						type: "loop",
						action: this,
						loopDelta: i
					})
				}
			} else this.time = n;
			if (a && 1 == (1 & r)) return t - n
		}
		return n
	}
	_setEndings(e, t, i) {
		const n = this._interpolantSettings;
		i ? (n.endingStart = 2401, n.endingEnd = 2401) : (n.endingStart = e ? this.zeroSlopeAtStart ? 2401 : 2400 : 2402, n.endingEnd = t ? this.zeroSlopeAtEnd ? 2401 : 2400 : 2402)
	}
	_scheduleFading(e, t, i) {
		const n = this._mixer,
			r = n.time;
		let a = this._weightInterpolant;
		null === a && (a = n._lendControlInterpolant(), this._weightInterpolant = a);
		const s = a.parameterPositions,
			o = a.sampleValues;
		return s[0] = r, o[0] = t, s[1] = r + e, o[1] = i, this
	}
}
class AnimationMixer extends EventDispatcher {
	constructor(e) {
		super(), this._root = e, this._initMemoryManager(), this._accuIndex = 0, this.time = 0, this.timeScale = 1
	}
	_bindAction(e, t) {
		const i = e._localRoot || this._root,
			n = e._clip.tracks,
			r = n.length,
			a = e._propertyBindings,
			s = e._interpolants,
			o = i.uuid,
			l = this._bindingsByRootAndName;
		let c = l[o];
		void 0 === c && (c = {}, l[o] = c);
		for (let e = 0; e !== r; ++e) {
			const r = n[e],
				l = r.name;
			let h = c[l];
			if (void 0 !== h) a[e] = h;
			else {
				if (h = a[e], void 0 !== h) {
					null === h._cacheIndex && (++h.referenceCount, this._addInactiveBinding(h, o, l));
					continue
				}
				const n = t && t._propertyBindings[e].binding.parsedPath;
				h = new PropertyMixer(PropertyBinding.create(i, l, n), r.ValueTypeName, r.getValueSize()), ++h.referenceCount, this._addInactiveBinding(h, o, l), a[e] = h
			}
			s[e].resultBuffer = h.buffer
		}
	}
	_activateAction(e) {
		if (!this._isActiveAction(e)) {
			if (null === e._cacheIndex) {
				const t = (e._localRoot || this._root).uuid,
					i = e._clip.uuid,
					n = this._actionsByClip[i];
				this._bindAction(e, n && n.knownActions[0]), this._addInactiveAction(e, i, t)
			}
			const t = e._propertyBindings;
			for (let e = 0, i = t.length; e !== i; ++e) {
				const i = t[e];
				0 == i.useCount++ && (this._lendBinding(i), i.saveOriginalState())
			}
			this._lendAction(e)
		}
	}
	_deactivateAction(e) {
		if (this._isActiveAction(e)) {
			const t = e._propertyBindings;
			for (let e = 0, i = t.length; e !== i; ++e) {
				const i = t[e];
				0 == --i.useCount && (i.restoreOriginalState(), this._takeBackBinding(i))
			}
			this._takeBackAction(e)
		}
	}
	_initMemoryManager() {
		this._actions = [], this._nActiveActions = 0, this._actionsByClip = {}, this._bindings = [], this._nActiveBindings = 0, this._bindingsByRootAndName = {}, this._controlInterpolants = [], this._nActiveControlInterpolants = 0;
		const e = this;
		this.stats = {
			actions: {
				get total() {
					return e._actions.length
				},
				get inUse() {
					return e._nActiveActions
				}
			},
			bindings: {
				get total() {
					return e._bindings.length
				},
				get inUse() {
					return e._nActiveBindings
				}
			},
			controlInterpolants: {
				get total() {
					return e._controlInterpolants.length
				},
				get inUse() {
					return e._nActiveControlInterpolants
				}
			}
		}
	}
	_isActiveAction(e) {
		const t = e._cacheIndex;
		return null !== t && t < this._nActiveActions
	}
	_addInactiveAction(e, t, i) {
		const n = this._actions,
			r = this._actionsByClip;
		let a = r[t];
		if (void 0 === a) a = {
			knownActions: [e],
			actionByRoot: {}
		}, e._byClipCacheIndex = 0, r[t] = a;
		else {
			const t = a.knownActions;
			e._byClipCacheIndex = t.length, t.push(e)
		}
		e._cacheIndex = n.length, n.push(e), a.actionByRoot[i] = e
	}
	_removeInactiveAction(e) {
		const t = this._actions,
			i = t[t.length - 1],
			n = e._cacheIndex;
		i._cacheIndex = n, t[n] = i, t.pop(), e._cacheIndex = null;
		const r = e._clip.uuid,
			a = this._actionsByClip,
			s = a[r],
			o = s.knownActions,
			l = o[o.length - 1],
			c = e._byClipCacheIndex;
		l._byClipCacheIndex = c, o[c] = l, o.pop(), e._byClipCacheIndex = null;
		delete s.actionByRoot[(e._localRoot || this._root).uuid], 0 === o.length && delete a[r], this._removeInactiveBindingsForAction(e)
	}
	_removeInactiveBindingsForAction(e) {
		const t = e._propertyBindings;
		for (let e = 0, i = t.length; e !== i; ++e) {
			const i = t[e];
			0 == --i.referenceCount && this._removeInactiveBinding(i)
		}
	}
	_lendAction(e) {
		const t = this._actions,
			i = e._cacheIndex,
			n = this._nActiveActions++,
			r = t[n];
		e._cacheIndex = n, t[n] = e, r._cacheIndex = i, t[i] = r
	}
	_takeBackAction(e) {
		const t = this._actions,
			i = e._cacheIndex,
			n = --this._nActiveActions,
			r = t[n];
		e._cacheIndex = n, t[n] = e, r._cacheIndex = i, t[i] = r
	}
	_addInactiveBinding(e, t, i) {
		const n = this._bindingsByRootAndName,
			r = this._bindings;
		let a = n[t];
		void 0 === a && (a = {}, n[t] = a), a[i] = e, e._cacheIndex = r.length, r.push(e)
	}
	_removeInactiveBinding(e) {
		const t = this._bindings,
			i = e.binding,
			n = i.rootNode.uuid,
			r = i.path,
			a = this._bindingsByRootAndName,
			s = a[n],
			o = t[t.length - 1],
			l = e._cacheIndex;
		o._cacheIndex = l, t[l] = o, t.pop(), delete s[r], 0 === Object.keys(s).length && delete a[n]
	}
	_lendBinding(e) {
		const t = this._bindings,
			i = e._cacheIndex,
			n = this._nActiveBindings++,
			r = t[n];
		e._cacheIndex = n, t[n] = e, r._cacheIndex = i, t[i] = r
	}
	_takeBackBinding(e) {
		const t = this._bindings,
			i = e._cacheIndex,
			n = --this._nActiveBindings,
			r = t[n];
		e._cacheIndex = n, t[n] = e, r._cacheIndex = i, t[i] = r
	}
	_lendControlInterpolant() {
		const e = this._controlInterpolants,
			t = this._nActiveControlInterpolants++;
		let i = e[t];
		return void 0 === i && (i = new LinearInterpolant(new Float32Array(2), new Float32Array(2), 1, this._controlInterpolantsResultBuffer), i.__cacheIndex = t, e[t] = i), i
	}
	_takeBackControlInterpolant(e) {
		const t = this._controlInterpolants,
			i = e.__cacheIndex,
			n = --this._nActiveControlInterpolants,
			r = t[n];
		e.__cacheIndex = n, t[n] = e, r.__cacheIndex = i, t[i] = r
	}
	clipAction(e, t, i) {
		const n = t || this._root,
			r = n.uuid;
		let a = "string" == typeof e ? AnimationClip.findByName(n, e) : e;
		const s = null !== a ? a.uuid : e,
			o = this._actionsByClip[s];
		let l = null;
		if (void 0 === i && (i = null !== a ? a.blendMode : 2500), void 0 !== o) {
			const e = o.actionByRoot[r];
			if (void 0 !== e && e.blendMode === i) return e;
			l = o.knownActions[0], null === a && (a = l._clip)
		}
		if (null === a) return null;
		const c = new AnimationAction(this, a, t, i);
		return this._bindAction(c, l), this._addInactiveAction(c, s, r), c
	}
	existingAction(e, t) {
		const i = t || this._root,
			n = i.uuid,
			r = "string" == typeof e ? AnimationClip.findByName(i, e) : e,
			a = r ? r.uuid : e,
			s = this._actionsByClip[a];
		return void 0 !== s && s.actionByRoot[n] || null
	}
	stopAllAction() {
		const e = this._actions;
		for (let t = this._nActiveActions - 1; t >= 0; --t) e[t].stop();
		return this
	}
	update(e) {
		e *= this.timeScale;
		const t = this._actions,
			i = this._nActiveActions,
			n = this.time += e,
			r = Math.sign(e),
			a = this._accuIndex ^= 1;
		for (let s = 0; s !== i; ++s) {
			t[s]._update(n, e, r, a)
		}
		const s = this._bindings,
			o = this._nActiveBindings;
		for (let e = 0; e !== o; ++e) s[e].apply(a);
		return this
	}
	setTime(e) {
		this.time = 0;
		for (let e = 0; e < this._actions.length; e++) this._actions[e].time = 0;
		return this.update(e)
	}
	getRoot() {
		return this._root
	}
	uncacheClip(e) {
		const t = this._actions,
			i = e.uuid,
			n = this._actionsByClip,
			r = n[i];
		if (void 0 !== r) {
			const e = r.knownActions;
			for (let i = 0, n = e.length; i !== n; ++i) {
				const n = e[i];
				this._deactivateAction(n);
				const r = n._cacheIndex,
					a = t[t.length - 1];
				n._cacheIndex = null, n._byClipCacheIndex = null, a._cacheIndex = r, t[r] = a, t.pop(), this._removeInactiveBindingsForAction(n)
			}
			delete n[i]
		}
	}
	uncacheRoot(e) {
		const t = e.uuid,
			i = this._actionsByClip;
		for (const e in i) {
			const n = i[e].actionByRoot[t];
			void 0 !== n && (this._deactivateAction(n), this._removeInactiveAction(n))
		}
		const n = this._bindingsByRootAndName[t];
		if (void 0 !== n)
			for (const e in n) {
				const t = n[e];
				t.restoreOriginalState(), this._removeInactiveBinding(t)
			}
	}
	uncacheAction(e, t) {
		const i = this.existingAction(e, t);
		null !== i && (this._deactivateAction(i), this._removeInactiveAction(i))
	}
}
AnimationMixer.prototype._controlInterpolantsResultBuffer = new Float32Array(1);
class Uniform {
	constructor(e) {
		"string" == typeof e && (console.warn("THREE.Uniform: Type parameter is no longer needed."), e = arguments[1]), this.value = e
	}
	clone() {
		return new Uniform(void 0 === this.value.clone ? this.value : this.value.clone())
	}
}

function InstancedInterleavedBuffer(e, t, i) {
	InterleavedBuffer.call(this, e, t), this.meshPerAttribute = i || 1
}

function GLBufferAttribute(e, t, i, n, r) {
	this.buffer = e, this.type = t, this.itemSize = i, this.elementSize = n, this.count = r, this.version = 0
}

function Raycaster(e, t, i = 0, n = 1 / 0) {
	this.ray = new Ray(e, t), this.near = i, this.far = n, this.camera = null, this.layers = new Layers, this.params = {
		Mesh: {},
		Line: {
			threshold: 1
		},
		LOD: {},
		Points: {
			threshold: 1
		},
		Sprite: {}
	}, Object.defineProperties(this.params, {
		PointCloud: {
			get: function() {
				return console.warn("THREE.Raycaster: params.PointCloud has been renamed to params.Points."), this.Points
			}
		}
	})
}

function ascSort(e, t) {
	return e.distance - t.distance
}

function intersectObject(e, t, i, n) {
	if (e.layers.test(t.layers) && e.raycast(t, i), !0 === n) {
		const n = e.children;
		for (let e = 0, r = n.length; e < r; e++) intersectObject(n[e], t, i, !0)
	}
}
InstancedInterleavedBuffer.prototype = Object.assign(Object.create(InterleavedBuffer.prototype), {
	constructor: InstancedInterleavedBuffer,
	isInstancedInterleavedBuffer: !0,
	copy: function(e) {
		return InterleavedBuffer.prototype.copy.call(this, e), this.meshPerAttribute = e.meshPerAttribute, this
	},
	clone: function(e) {
		const t = InterleavedBuffer.prototype.clone.call(this, e);
		return t.meshPerAttribute = this.meshPerAttribute, t
	},
	toJSON: function(e) {
		const t = InterleavedBuffer.prototype.toJSON.call(this, e);
		return t.isInstancedInterleavedBuffer = !0, t.meshPerAttribute = this.meshPerAttribute, t
	}
}), Object.defineProperty(GLBufferAttribute.prototype, "needsUpdate", {
	set: function(e) {
		!0 === e && this.version++
	}
}), Object.assign(GLBufferAttribute.prototype, {
	isGLBufferAttribute: !0,
	setBuffer: function(e) {
		return this.buffer = e, this
	},
	setType: function(e, t) {
		return this.type = e, this.elementSize = t, this
	},
	setItemSize: function(e) {
		return this.itemSize = e, this
	},
	setCount: function(e) {
		return this.count = e, this
	}
}), Object.assign(Raycaster.prototype, {
	set: function(e, t) {
		this.ray.set(e, t)
	},
	setFromCamera: function(e, t) {
		t && t.isPerspectiveCamera ? (this.ray.origin.setFromMatrixPosition(t.matrixWorld), this.ray.direction.set(e.x, e.y, .5).unproject(t).sub(this.ray.origin).normalize(), this.camera = t) : t && t.isOrthographicCamera ? (this.ray.origin.set(e.x, e.y, (t.near + t.far) / (t.near - t.far)).unproject(t), this.ray.direction.set(0, 0, -1).transformDirection(t.matrixWorld), this.camera = t) : console.error("THREE.Raycaster: Unsupported camera type: " + t.type)
	},
	intersectObject: function(e, t = !1, i = []) {
		return intersectObject(e, this, i, t), i.sort(ascSort), i
	},
	intersectObjects: function(e, t = !1, i = []) {
		for (let n = 0, r = e.length; n < r; n++) intersectObject(e[n], this, i, t);
		return i.sort(ascSort), i
	}
});
class Spherical {
	constructor(e = 1, t = 0, i = 0) {
		return this.radius = e, this.phi = t, this.theta = i, this
	}
	set(e, t, i) {
		return this.radius = e, this.phi = t, this.theta = i, this
	}
	copy(e) {
		return this.radius = e.radius, this.phi = e.phi, this.theta = e.theta, this
	}
	makeSafe() {
		return this.phi = Math.max(1e-6, Math.min(Math.PI - 1e-6, this.phi)), this
	}
	setFromVector3(e) {
		return this.setFromCartesianCoords(e.x, e.y, e.z)
	}
	setFromCartesianCoords(e, t, i) {
		return this.radius = Math.sqrt(e * e + t * t + i * i), 0 === this.radius ? (this.theta = 0, this.phi = 0) : (this.theta = Math.atan2(e, i), this.phi = Math.acos(MathUtils.clamp(t / this.radius, -1, 1))), this
	}
	clone() {
		return (new this.constructor).copy(this)
	}
}

function ImmediateRenderObject(e) {
	Object3D.call(this), this.material = e, this.render = function() {}, this.hasPositions = !1, this.hasNormals = !1, this.hasColors = !1, this.hasUvs = !1, this.positionArray = null, this.normalArray = null, this.colorArray = null, this.uvArray = null, this.count = 0
}
ImmediateRenderObject.prototype = Object.create(Object3D.prototype), ImmediateRenderObject.prototype.constructor = ImmediateRenderObject, ImmediateRenderObject.prototype.isImmediateRenderObject = !0;
const _vector$a = new Vector3,
	_boneMatrix = new Matrix4,
	_matrixWorldInv = new Matrix4;
class SkeletonHelper extends LineSegments {
	constructor(e) {
		const t = getBoneList(e),
			i = new BufferGeometry,
			n = [],
			r = [],
			a = new Color(0, 0, 1),
			s = new Color(0, 1, 0);
		for (let e = 0; e < t.length; e++) {
			const i = t[e];
			i.parent && i.parent.isBone && (n.push(0, 0, 0), n.push(0, 0, 0), r.push(a.r, a.g, a.b), r.push(s.r, s.g, s.b))
		}
		i.setAttribute("position", new Float32BufferAttribute(n, 3)), i.setAttribute("color", new Float32BufferAttribute(r, 3));
		super(i, new LineBasicMaterial({
			vertexColors: !0,
			depthTest: !1,
			depthWrite: !1,
			toneMapped: !1,
			transparent: !0
		})), this.type = "SkeletonHelper", this.isSkeletonHelper = !0, this.root = e, this.bones = t, this.matrix = e.matrixWorld, this.matrixAutoUpdate = !1
	}
	updateMatrixWorld(e) {
		const t = this.bones,
			i = this.geometry,
			n = i.getAttribute("position");
		_matrixWorldInv.copy(this.root.matrixWorld).invert();
		for (let e = 0, i = 0; e < t.length; e++) {
			const r = t[e];
			r.parent && r.parent.isBone && (_boneMatrix.multiplyMatrices(_matrixWorldInv, r.matrixWorld), _vector$a.setFromMatrixPosition(_boneMatrix), n.setXYZ(i, _vector$a.x, _vector$a.y, _vector$a.z), _boneMatrix.multiplyMatrices(_matrixWorldInv, r.parent.matrixWorld), _vector$a.setFromMatrixPosition(_boneMatrix), n.setXYZ(i + 1, _vector$a.x, _vector$a.y, _vector$a.z), i += 2)
		}
		i.getAttribute("position").needsUpdate = !0, super.updateMatrixWorld(e)
	}
}

function getBoneList(e) {
	const t = [];
	e && e.isBone && t.push(e);
	for (let i = 0; i < e.children.length; i++) t.push.apply(t, getBoneList(e.children[i]));
	return t
}
const _floatView = new Float32Array(1),
	_int32View = new Int32Array(_floatView.buffer),
	DataUtils = {
		toHalfFloat: function(e) {
			_floatView[0] = e;
			const t = _int32View[0];
			let i = t >> 16 & 32768,
				n = t >> 12 & 2047;
			const r = t >> 23 & 255;
			return r < 103 ? i : r > 142 ? (i |= 31744, i |= (255 == r ? 0 : 1) && 8388607 & t, i) : r < 113 ? (n |= 2048, i |= (n >> 114 - r) + (n >> 113 - r & 1), i) : (i |= r - 112 << 10 | n >> 1, i += 1 & n, i)
		}
	},
	LOD_MIN = 4,
	LOD_MAX = 8,
	SIZE_MAX = Math.pow(2, 8),
	EXTRA_LOD_SIGMA = [.125, .215, .35, .446, .526, .582],
	TOTAL_LODS = 5 + EXTRA_LOD_SIGMA.length,
	MAX_SAMPLES = 20,
	ENCODINGS = {
		3e3: 0,
		3001: 1,
		3002: 2,
		3004: 3,
		3005: 4,
		3006: 5,
		3007: 6
	},
	backgroundMaterial = new MeshBasicMaterial({
		side: 1,
		depthWrite: !1,
		depthTest: !1
	}),
	backgroundBox = new Mesh(new BoxGeometry, backgroundMaterial),
	_flatCamera$1 = new OrthographicCamera,
	{
		_lodPlanes: _lodPlanes,
		_sizeLods: _sizeLods,
		_sigmas: _sigmas
	} = _createPlanes(),
	_clearColor = new Color;
let _oldTarget = null;
const PHI = (1 + Math.sqrt(5)) / 2,
	INV_PHI = 1 / PHI,
	_axisDirections = [new Vector3(1, 1, 1), new Vector3(-1, 1, 1), new Vector3(1, 1, -1), new Vector3(-1, 1, -1), new Vector3(0, PHI, INV_PHI), new Vector3(0, PHI, -INV_PHI), new Vector3(INV_PHI, 0, PHI), new Vector3(-INV_PHI, 0, PHI), new Vector3(PHI, INV_PHI, 0), new Vector3(-PHI, INV_PHI, 0)];

function convertLinearToRGBE(e) {
	const t = Math.max(e.r, e.g, e.b),
		i = Math.min(Math.max(Math.ceil(Math.log2(t)), -128), 127);
	e.multiplyScalar(Math.pow(2, -i));
	return (i + 128) / 255
}
class PMREMGenerator {
	constructor(e) {
		this._renderer = e, this._pingPongRenderTarget = null, this._blurMaterial = _getBlurShader(20), this._equirectShader = null, this._cubemapShader = null, this._compileMaterial(this._blurMaterial)
	}
	fromScene(e, t = 0, i = .1, n = 100) {
		_oldTarget = this._renderer.getRenderTarget();
		const r = this._allocateTargets();
		return this._sceneToCubeUV(e, i, n, r), t > 0 && this._blur(r, 0, 0, t), this._applyPMREM(r), this._cleanup(r), r
	}
	fromEquirectangular(e) {
		return this._fromTexture(e)
	}
	fromCubemap(e) {
		return this._fromTexture(e)
	}
	compileCubemapShader() {
		null === this._cubemapShader && (this._cubemapShader = _getCubemapShader(), this._compileMaterial(this._cubemapShader))
	}
	compileEquirectangularShader() {
		null === this._equirectShader && (this._equirectShader = _getEquirectShader(), this._compileMaterial(this._equirectShader))
	}
	dispose() {
		this._blurMaterial.dispose(), null !== this._cubemapShader && this._cubemapShader.dispose(), null !== this._equirectShader && this._equirectShader.dispose();
		for (let e = 0; e < _lodPlanes.length; e++) _lodPlanes[e].dispose()
	}
	_cleanup(e) {
		this._pingPongRenderTarget.dispose(), this._renderer.setRenderTarget(_oldTarget), e.scissorTest = !1, _setViewport(e, 0, 0, e.width, e.height)
	}
	_fromTexture(e) {
		_oldTarget = this._renderer.getRenderTarget();
		const t = this._allocateTargets(e);
		return this._textureToCubeUV(e, t), this._applyPMREM(t), this._cleanup(t), t
	}
	_allocateTargets(e) {
		const t = {
				magFilter: 1003,
				minFilter: 1003,
				generateMipmaps: !1,
				type: 1009,
				format: 1023,
				encoding: _isLDR(e) ? e.encoding : 3002,
				depthBuffer: !1
			},
			i = _createRenderTarget(t);
		return i.depthBuffer = !e, this._pingPongRenderTarget = _createRenderTarget(t), i
	}
	_compileMaterial(e) {
		const t = new Mesh(_lodPlanes[0], e);
		this._renderer.compile(t, _flatCamera$1)
	}
	_sceneToCubeUV(e, t, i, n) {
		const r = new PerspectiveCamera(90, 1, t, i),
			a = [1, -1, 1, 1, 1, 1],
			s = [1, 1, 1, -1, -1, -1],
			o = this._renderer,
			l = o.autoClear,
			c = o.outputEncoding,
			h = o.toneMapping;
		o.getClearColor(_clearColor), o.toneMapping = 0, o.outputEncoding = 3e3, o.autoClear = !1;
		let u = !1;
		const d = e.background;
		if (d) {
			if (d.isColor) {
				backgroundMaterial.color.copy(d).convertSRGBToLinear(), e.background = null;
				const t = convertLinearToRGBE(backgroundMaterial.color);
				backgroundMaterial.opacity = t, u = !0
			}
		} else {
			backgroundMaterial.color.copy(_clearColor).convertSRGBToLinear();
			const e = convertLinearToRGBE(backgroundMaterial.color);
			backgroundMaterial.opacity = e, u = !0
		}
		for (let t = 0; t < 6; t++) {
			const i = t % 3;
			0 == i ? (r.up.set(0, a[t], 0), r.lookAt(s[t], 0, 0)) : 1 == i ? (r.up.set(0, 0, a[t]), r.lookAt(0, s[t], 0)) : (r.up.set(0, a[t], 0), r.lookAt(0, 0, s[t])), _setViewport(n, i * SIZE_MAX, t > 2 ? SIZE_MAX : 0, SIZE_MAX, SIZE_MAX), o.setRenderTarget(n), u && o.render(backgroundBox, r), o.render(e, r)
		}
		o.toneMapping = h, o.outputEncoding = c, o.autoClear = l
	}
	_textureToCubeUV(e, t) {
		const i = this._renderer;
		e.isCubeTexture ? null == this._cubemapShader && (this._cubemapShader = _getCubemapShader()) : null == this._equirectShader && (this._equirectShader = _getEquirectShader());
		const n = e.isCubeTexture ? this._cubemapShader : this._equirectShader,
			r = new Mesh(_lodPlanes[0], n),
			a = n.uniforms;
		a.envMap.value = e, e.isCubeTexture || a.texelSize.value.set(1 / e.image.width, 1 / e.image.height), a.inputEncoding.value = ENCODINGS[e.encoding], a.outputEncoding.value = ENCODINGS[t.texture.encoding], _setViewport(t, 0, 0, 3 * SIZE_MAX, 2 * SIZE_MAX), i.setRenderTarget(t), i.render(r, _flatCamera$1)
	}
	_applyPMREM(e) {
		const t = this._renderer,
			i = t.autoClear;
		t.autoClear = !1;
		for (let t = 1; t < TOTAL_LODS; t++) {
			const i = Math.sqrt(_sigmas[t] * _sigmas[t] - _sigmas[t - 1] * _sigmas[t - 1]),
				n = _axisDirections[(t - 1) % _axisDirections.length];
			this._blur(e, t - 1, t, i, n)
		}
		t.autoClear = i
	}
	_blur(e, t, i, n, r) {
		const a = this._pingPongRenderTarget;
		this._halfBlur(e, a, t, i, n, "latitudinal", r), this._halfBlur(a, e, i, i, n, "longitudinal", r)
	}
	_halfBlur(e, t, i, n, r, a, s) {
		const o = this._renderer,
			l = this._blurMaterial;
		"latitudinal" !== a && "longitudinal" !== a && console.error("blur direction must be either latitudinal or longitudinal!");
		const c = new Mesh(_lodPlanes[n], l),
			h = l.uniforms,
			u = _sizeLods[i] - 1,
			d = isFinite(r) ? Math.PI / (2 * u) : 2 * Math.PI / 39,
			p = r / d,
			m = isFinite(r) ? 1 + Math.floor(3 * p) : 20;
		m > 20 && console.warn(`sigmaRadians, ${r}, is too large and will clip, as it requested ${m} samples when the maximum is set to 20`);
		const A = [];
		let g = 0;
		for (let e = 0; e < 20; ++e) {
			const t = e / p,
				i = Math.exp(-t * t / 2);
			A.push(i), 0 == e ? g += i : e < m && (g += 2 * i)
		}
		for (let e = 0; e < A.length; e++) A[e] = A[e] / g;
		h.envMap.value = e.texture, h.samples.value = m, h.weights.value = A, h.latitudinal.value = "latitudinal" === a, s && (h.poleAxis.value = s), h.dTheta.value = d, h.mipInt.value = 8 - i, h.inputEncoding.value = ENCODINGS[e.texture.encoding], h.outputEncoding.value = ENCODINGS[e.texture.encoding];
		const f = _sizeLods[n];
		_setViewport(t, 3 * Math.max(0, SIZE_MAX - 2 * f), (0 === n ? 0 : 2 * SIZE_MAX) + 2 * f * (n > 4 ? n - 8 + 4 : 0), 3 * f, 2 * f), o.setRenderTarget(t), o.render(c, _flatCamera$1)
	}
}

function _isLDR(e) {
	return void 0 !== e && 1009 === e.type && (3e3 === e.encoding || 3001 === e.encoding || 3007 === e.encoding)
}

function _createPlanes() {
	const e = [],
		t = [],
		i = [];
	let n = 8;
	for (let r = 0; r < TOTAL_LODS; r++) {
		const a = Math.pow(2, n);
		t.push(a);
		let s = 1 / a;
		r > 4 ? s = EXTRA_LOD_SIGMA[r - 8 + 4 - 1] : 0 == r && (s = 0), i.push(s);
		const o = 1 / (a - 1),
			l = -o / 2,
			c = 1 + o / 2,
			h = [l, l, c, l, c, c, l, l, c, c, l, c],
			u = 6,
			d = 6,
			p = 3,
			m = 2,
			A = 1,
			g = new Float32Array(p * d * u),
			f = new Float32Array(m * d * u),
			v = new Float32Array(A * d * u);
		for (let e = 0; e < u; e++) {
			const t = e % 3 * 2 / 3 - 1,
				i = e > 2 ? 0 : -1,
				n = [t, i, 0, t + 2 / 3, i, 0, t + 2 / 3, i + 1, 0, t, i, 0, t + 2 / 3, i + 1, 0, t, i + 1, 0];
			g.set(n, p * d * e), f.set(h, m * d * e);
			const r = [e, e, e, e, e, e];
			v.set(r, A * d * e)
		}
		const y = new BufferGeometry;
		y.setAttribute("position", new BufferAttribute(g, p)), y.setAttribute("uv", new BufferAttribute(f, m)), y.setAttribute("faceIndex", new BufferAttribute(v, A)), e.push(y), n > 4 && n--
	}
	return {
		_lodPlanes: e,
		_sizeLods: t,
		_sigmas: i
	}
}

function _createRenderTarget(e) {
	const t = new WebGLRenderTarget(3 * SIZE_MAX, 3 * SIZE_MAX, e);
	return t.texture.mapping = 306, t.texture.name = "PMREM.cubeUv", t.scissorTest = !0, t
}

function _setViewport(e, t, i, n, r) {
	e.viewport.set(t, i, n, r), e.scissor.set(t, i, n, r)
}

function _getBlurShader(e) {
	const t = new Float32Array(e),
		i = new Vector3(0, 1, 0);
	return new RawShaderMaterial({
		name: "SphericalGaussianBlur",
		defines: {
			n: e
		},
		uniforms: {
			envMap: {
				value: null
			},
			samples: {
				value: 1
			},
			weights: {
				value: t
			},
			latitudinal: {
				value: !1
			},
			dTheta: {
				value: 0
			},
			mipInt: {
				value: 0
			},
			poleAxis: {
				value: i
			},
			inputEncoding: {
				value: ENCODINGS[3e3]
			},
			outputEncoding: {
				value: ENCODINGS[3e3]
			}
		},
		vertexShader: _getCommonVertexShader(),
		fragmentShader: `\n\n\t\t\tprecision mediump float;\n\t\t\tprecision mediump int;\n\n\t\t\tvarying vec3 vOutputDirection;\n\n\t\t\tuniform sampler2D envMap;\n\t\t\tuniform int samples;\n\t\t\tuniform float weights[ n ];\n\t\t\tuniform bool latitudinal;\n\t\t\tuniform float dTheta;\n\t\t\tuniform float mipInt;\n\t\t\tuniform vec3 poleAxis;\n\n\t\t\t${_getEncodings()}\n\n\t\t\t#define ENVMAP_TYPE_CUBE_UV\n\t\t\t#include <cube_uv_reflection_fragment>\n\n\t\t\tvec3 getSample( float theta, vec3 axis ) {\n\n\t\t\t\tfloat cosTheta = cos( theta );\n\t\t\t\t// Rodrigues' axis-angle rotation\n\t\t\t\tvec3 sampleDirection = vOutputDirection * cosTheta\n\t\t\t\t\t+ cross( axis, vOutputDirection ) * sin( theta )\n\t\t\t\t\t+ axis * dot( axis, vOutputDirection ) * ( 1.0 - cosTheta );\n\n\t\t\t\treturn bilinearCubeUV( envMap, sampleDirection, mipInt );\n\n\t\t\t}\n\n\t\t\tvoid main() {\n\n\t\t\t\tvec3 axis = latitudinal ? poleAxis : cross( poleAxis, vOutputDirection );\n\n\t\t\t\tif ( all( equal( axis, vec3( 0.0 ) ) ) ) {\n\n\t\t\t\t\taxis = vec3( vOutputDirection.z, 0.0, - vOutputDirection.x );\n\n\t\t\t\t}\n\n\t\t\t\taxis = normalize( axis );\n\n\t\t\t\tgl_FragColor = vec4( 0.0, 0.0, 0.0, 1.0 );\n\t\t\t\tgl_FragColor.rgb += weights[ 0 ] * getSample( 0.0, axis );\n\n\t\t\t\tfor ( int i = 1; i < n; i++ ) {\n\n\t\t\t\t\tif ( i >= samples ) {\n\n\t\t\t\t\t\tbreak;\n\n\t\t\t\t\t}\n\n\t\t\t\t\tfloat theta = dTheta * float( i );\n\t\t\t\t\tgl_FragColor.rgb += weights[ i ] * getSample( -1.0 * theta, axis );\n\t\t\t\t\tgl_FragColor.rgb += weights[ i ] * getSample( theta, axis );\n\n\t\t\t\t}\n\n\t\t\t\tgl_FragColor = linearToOutputTexel( gl_FragColor );\n\n\t\t\t}\n\t\t`,
		blending: 0,
		depthTest: !1,
		depthWrite: !1
	})
}

function _getEquirectShader() {
	const e = new Vector2(1, 1);
	return new RawShaderMaterial({
		name: "EquirectangularToCubeUV",
		uniforms: {
			envMap: {
				value: null
			},
			texelSize: {
				value: e
			},
			inputEncoding: {
				value: ENCODINGS[3e3]
			},
			outputEncoding: {
				value: ENCODINGS[3e3]
			}
		},
		vertexShader: _getCommonVertexShader(),
		fragmentShader: `\n\n\t\t\tprecision mediump float;\n\t\t\tprecision mediump int;\n\n\t\t\tvarying vec3 vOutputDirection;\n\n\t\t\tuniform sampler2D envMap;\n\t\t\tuniform vec2 texelSize;\n\n\t\t\t${_getEncodings()}\n\n\t\t\t#include <common>\n\n\t\t\tvoid main() {\n\n\t\t\t\tgl_FragColor = vec4( 0.0, 0.0, 0.0, 1.0 );\n\n\t\t\t\tvec3 outputDirection = normalize( vOutputDirection );\n\t\t\t\tvec2 uv = equirectUv( outputDirection );\n\n\t\t\t\tvec2 f = fract( uv / texelSize - 0.5 );\n\t\t\t\tuv -= f * texelSize;\n\t\t\t\tvec3 tl = envMapTexelToLinear( texture2D ( envMap, uv ) ).rgb;\n\t\t\t\tuv.x += texelSize.x;\n\t\t\t\tvec3 tr = envMapTexelToLinear( texture2D ( envMap, uv ) ).rgb;\n\t\t\t\tuv.y += texelSize.y;\n\t\t\t\tvec3 br = envMapTexelToLinear( texture2D ( envMap, uv ) ).rgb;\n\t\t\t\tuv.x -= texelSize.x;\n\t\t\t\tvec3 bl = envMapTexelToLinear( texture2D ( envMap, uv ) ).rgb;\n\n\t\t\t\tvec3 tm = mix( tl, tr, f.x );\n\t\t\t\tvec3 bm = mix( bl, br, f.x );\n\t\t\t\tgl_FragColor.rgb = mix( tm, bm, f.y );\n\n\t\t\t\tgl_FragColor = linearToOutputTexel( gl_FragColor );\n\n\t\t\t}\n\t\t`,
		blending: 0,
		depthTest: !1,
		depthWrite: !1
	})
}

function _getCubemapShader() {
	return new RawShaderMaterial({
		name: "CubemapToCubeUV",
		uniforms: {
			envMap: {
				value: null
			},
			inputEncoding: {
				value: ENCODINGS[3e3]
			},
			outputEncoding: {
				value: ENCODINGS[3e3]
			}
		},
		vertexShader: _getCommonVertexShader(),
		fragmentShader: `\n\n\t\t\tprecision mediump float;\n\t\t\tprecision mediump int;\n\n\t\t\tvarying vec3 vOutputDirection;\n\n\t\t\tuniform samplerCube envMap;\n\n\t\t\t${_getEncodings()}\n\n\t\t\tvoid main() {\n\n\t\t\t\tgl_FragColor = vec4( 0.0, 0.0, 0.0, 1.0 );\n\t\t\t\tgl_FragColor.rgb = envMapTexelToLinear( textureCube( envMap, vec3( - vOutputDirection.x, vOutputDirection.yz ) ) ).rgb;\n\t\t\t\tgl_FragColor = linearToOutputTexel( gl_FragColor );\n\n\t\t\t}\n\t\t`,
		blending: 0,
		depthTest: !1,
		depthWrite: !1
	})
}

function _getCommonVertexShader() {
	return "\n\n\t\tprecision mediump float;\n\t\tprecision mediump int;\n\n\t\tattribute vec3 position;\n\t\tattribute vec2 uv;\n\t\tattribute float faceIndex;\n\n\t\tvarying vec3 vOutputDirection;\n\n\t\t// RH coordinate system; PMREM face-indexing convention\n\t\tvec3 getDirection( vec2 uv, float face ) {\n\n\t\t\tuv = 2.0 * uv - 1.0;\n\n\t\t\tvec3 direction = vec3( uv, 1.0 );\n\n\t\t\tif ( face == 0.0 ) {\n\n\t\t\t\tdirection = direction.zyx; // ( 1, v, u ) pos x\n\n\t\t\t} else if ( face == 1.0 ) {\n\n\t\t\t\tdirection = direction.xzy;\n\t\t\t\tdirection.xz *= -1.0; // ( -u, 1, -v ) pos y\n\n\t\t\t} else if ( face == 2.0 ) {\n\n\t\t\t\tdirection.x *= -1.0; // ( -u, v, 1 ) pos z\n\n\t\t\t} else if ( face == 3.0 ) {\n\n\t\t\t\tdirection = direction.zyx;\n\t\t\t\tdirection.xz *= -1.0; // ( -1, v, -u ) neg x\n\n\t\t\t} else if ( face == 4.0 ) {\n\n\t\t\t\tdirection = direction.xzy;\n\t\t\t\tdirection.xy *= -1.0; // ( -u, -1, v ) neg y\n\n\t\t\t} else if ( face == 5.0 ) {\n\n\t\t\t\tdirection.z *= -1.0; // ( u, v, -1 ) neg z\n\n\t\t\t}\n\n\t\t\treturn direction;\n\n\t\t}\n\n\t\tvoid main() {\n\n\t\t\tvOutputDirection = getDirection( uv, faceIndex );\n\t\t\tgl_Position = vec4( position, 1.0 );\n\n\t\t}\n\t"
}

function _getEncodings() {
	return "\n\n\t\tuniform int inputEncoding;\n\t\tuniform int outputEncoding;\n\n\t\t#include <encodings_pars_fragment>\n\n\t\tvec4 inputTexelToLinear( vec4 value ) {\n\n\t\t\tif ( inputEncoding == 0 ) {\n\n\t\t\t\treturn value;\n\n\t\t\t} else if ( inputEncoding == 1 ) {\n\n\t\t\t\treturn sRGBToLinear( value );\n\n\t\t\t} else if ( inputEncoding == 2 ) {\n\n\t\t\t\treturn RGBEToLinear( value );\n\n\t\t\t} else if ( inputEncoding == 3 ) {\n\n\t\t\t\treturn RGBMToLinear( value, 7.0 );\n\n\t\t\t} else if ( inputEncoding == 4 ) {\n\n\t\t\t\treturn RGBMToLinear( value, 16.0 );\n\n\t\t\t} else if ( inputEncoding == 5 ) {\n\n\t\t\t\treturn RGBDToLinear( value, 256.0 );\n\n\t\t\t} else {\n\n\t\t\t\treturn GammaToLinear( value, 2.2 );\n\n\t\t\t}\n\n\t\t}\n\n\t\tvec4 linearToOutputTexel( vec4 value ) {\n\n\t\t\tif ( outputEncoding == 0 ) {\n\n\t\t\t\treturn value;\n\n\t\t\t} else if ( outputEncoding == 1 ) {\n\n\t\t\t\treturn LinearTosRGB( value );\n\n\t\t\t} else if ( outputEncoding == 2 ) {\n\n\t\t\t\treturn LinearToRGBE( value );\n\n\t\t\t} else if ( outputEncoding == 3 ) {\n\n\t\t\t\treturn LinearToRGBM( value, 7.0 );\n\n\t\t\t} else if ( outputEncoding == 4 ) {\n\n\t\t\t\treturn LinearToRGBM( value, 16.0 );\n\n\t\t\t} else if ( outputEncoding == 5 ) {\n\n\t\t\t\treturn LinearToRGBD( value, 256.0 );\n\n\t\t\t} else {\n\n\t\t\t\treturn LinearToGamma( value, 2.2 );\n\n\t\t\t}\n\n\t\t}\n\n\t\tvec4 envMapTexelToLinear( vec4 color ) {\n\n\t\t\treturn inputTexelToLinear( color );\n\n\t\t}\n\t"
}
Curve.create = function(e, t) {
	return console.log("THREE.Curve.create() has been deprecated"), e.prototype = Object.create(Curve.prototype), e.prototype.constructor = e, e.prototype.getPoint = t, e
}, SkeletonHelper.prototype.update = function() {
	console.error("THREE.SkeletonHelper: update() no longer needs to be called.")
}, Loader.prototype.extractUrlBase = function(e) {
	return console.warn("THREE.Loader: .extractUrlBase() has been deprecated. Use THREE.LoaderUtils.extractUrlBase() instead."), LoaderUtils.extractUrlBase(e)
}, Loader.Handlers = {
	add: function() {
		console.error("THREE.Loader: Handlers.add() has been removed. Use LoadingManager.addHandler() instead.")
	},
	get: function() {
		console.error("THREE.Loader: Handlers.get() has been removed. Use LoadingManager.getHandler() instead.")
	}
}, Box3.prototype.center = function(e) {
	return console.warn("THREE.Box3: .center() has been renamed to .getCenter()."), this.getCenter(e)
}, Box3.prototype.empty = function() {
	return console.warn("THREE.Box3: .empty() has been renamed to .isEmpty()."), this.isEmpty()
}, Box3.prototype.isIntersectionBox = function(e) {
	return console.warn("THREE.Box3: .isIntersectionBox() has been renamed to .intersectsBox()."), this.intersectsBox(e)
}, Box3.prototype.isIntersectionSphere = function(e) {
	return console.warn("THREE.Box3: .isIntersectionSphere() has been renamed to .intersectsSphere()."), this.intersectsSphere(e)
}, Box3.prototype.size = function(e) {
	return console.warn("THREE.Box3: .size() has been renamed to .getSize()."), this.getSize(e)
}, Sphere.prototype.empty = function() {
	return console.warn("THREE.Sphere: .empty() has been renamed to .isEmpty()."), this.isEmpty()
}, Frustum.prototype.setFromMatrix = function(e) {
	return console.warn("THREE.Frustum: .setFromMatrix() has been renamed to .setFromProjectionMatrix()."), this.setFromProjectionMatrix(e)
}, MathUtils.random16 = function() {
	return console.warn("THREE.Math: .random16() has been deprecated. Use Math.random() instead."), Math.random()
}, MathUtils.nearestPowerOfTwo = function(e) {
	return console.warn("THREE.Math: .nearestPowerOfTwo() has been renamed to .floorPowerOfTwo()."), MathUtils.floorPowerOfTwo(e)
}, MathUtils.nextPowerOfTwo = function(e) {
	return console.warn("THREE.Math: .nextPowerOfTwo() has been renamed to .ceilPowerOfTwo()."), MathUtils.ceilPowerOfTwo(e)
}, Matrix3.prototype.flattenToArrayOffset = function(e, t) {
	return console.warn("THREE.Matrix3: .flattenToArrayOffset() has been deprecated. Use .toArray() instead."), this.toArray(e, t)
}, Matrix3.prototype.multiplyVector3 = function(e) {
	return console.warn("THREE.Matrix3: .multiplyVector3() has been removed. Use vector.applyMatrix3( matrix ) instead."), e.applyMatrix3(this)
}, Matrix3.prototype.multiplyVector3Array = function() {
	console.error("THREE.Matrix3: .multiplyVector3Array() has been removed.")
}, Matrix3.prototype.applyToBufferAttribute = function(e) {
	return console.warn("THREE.Matrix3: .applyToBufferAttribute() has been removed. Use attribute.applyMatrix3( matrix ) instead."), e.applyMatrix3(this)
}, Matrix3.prototype.applyToVector3Array = function() {
	console.error("THREE.Matrix3: .applyToVector3Array() has been removed.")
}, Matrix3.prototype.getInverse = function(e) {
	return console.warn("THREE.Matrix3: .getInverse() has been removed. Use matrixInv.copy( matrix ).invert(); instead."), this.copy(e).invert()
}, Matrix4.prototype.extractPosition = function(e) {
	return console.warn("THREE.Matrix4: .extractPosition() has been renamed to .copyPosition()."), this.copyPosition(e)
}, Matrix4.prototype.flattenToArrayOffset = function(e, t) {
	return console.warn("THREE.Matrix4: .flattenToArrayOffset() has been deprecated. Use .toArray() instead."), this.toArray(e, t)
}, Matrix4.prototype.getPosition = function() {
	return console.warn("THREE.Matrix4: .getPosition() has been removed. Use Vector3.setFromMatrixPosition( matrix ) instead."), (new Vector3).setFromMatrixColumn(this, 3)
}, Matrix4.prototype.setRotationFromQuaternion = function(e) {
	return console.warn("THREE.Matrix4: .setRotationFromQuaternion() has been renamed to .makeRotationFromQuaternion()."), this.makeRotationFromQuaternion(e)
}, Matrix4.prototype.multiplyToArray = function() {
	console.warn("THREE.Matrix4: .multiplyToArray() has been removed.")
}, Matrix4.prototype.multiplyVector3 = function(e) {
	return console.warn("THREE.Matrix4: .multiplyVector3() has been removed. Use vector.applyMatrix4( matrix ) instead."), e.applyMatrix4(this)
}, Matrix4.prototype.multiplyVector4 = function(e) {
	return console.warn("THREE.Matrix4: .multiplyVector4() has been removed. Use vector.applyMatrix4( matrix ) instead."), e.applyMatrix4(this)
}, Matrix4.prototype.multiplyVector3Array = function() {
	console.error("THREE.Matrix4: .multiplyVector3Array() has been removed.")
}, Matrix4.prototype.rotateAxis = function(e) {
	console.warn("THREE.Matrix4: .rotateAxis() has been removed. Use Vector3.transformDirection( matrix ) instead."), e.transformDirection(this)
}, Matrix4.prototype.crossVector = function(e) {
	return console.warn("THREE.Matrix4: .crossVector() has been removed. Use vector.applyMatrix4( matrix ) instead."), e.applyMatrix4(this)
}, Matrix4.prototype.translate = function() {
	console.error("THREE.Matrix4: .translate() has been removed.")
}, Matrix4.prototype.rotateX = function() {
	console.error("THREE.Matrix4: .rotateX() has been removed.")
}, Matrix4.prototype.rotateY = function() {
	console.error("THREE.Matrix4: .rotateY() has been removed.")
}, Matrix4.prototype.rotateZ = function() {
	console.error("THREE.Matrix4: .rotateZ() has been removed.")
}, Matrix4.prototype.rotateByAxis = function() {
	console.error("THREE.Matrix4: .rotateByAxis() has been removed.")
}, Matrix4.prototype.applyToBufferAttribute = function(e) {
	return console.warn("THREE.Matrix4: .applyToBufferAttribute() has been removed. Use attribute.applyMatrix4( matrix ) instead."), e.applyMatrix4(this)
}, Matrix4.prototype.applyToVector3Array = function() {
	console.error("THREE.Matrix4: .applyToVector3Array() has been removed.")
}, Matrix4.prototype.makeFrustum = function(e, t, i, n, r, a) {
	return console.warn("THREE.Matrix4: .makeFrustum() has been removed. Use .makePerspective( left, right, top, bottom, near, far ) instead."), this.makePerspective(e, t, n, i, r, a)
}, Matrix4.prototype.getInverse = function(e) {
	return console.warn("THREE.Matrix4: .getInverse() has been removed. Use matrixInv.copy( matrix ).invert(); instead."), this.copy(e).invert()
}, Plane.prototype.isIntersectionLine = function(e) {
	return console.warn("THREE.Plane: .isIntersectionLine() has been renamed to .intersectsLine()."), this.intersectsLine(e)
}, Quaternion.prototype.multiplyVector3 = function(e) {
	return console.warn("THREE.Quaternion: .multiplyVector3() has been removed. Use is now vector.applyQuaternion( quaternion ) instead."), e.applyQuaternion(this)
}, Quaternion.prototype.inverse = function() {
	return console.warn("THREE.Quaternion: .inverse() has been renamed to invert()."), this.invert()
}, Ray.prototype.isIntersectionBox = function(e) {
	return console.warn("THREE.Ray: .isIntersectionBox() has been renamed to .intersectsBox()."), this.intersectsBox(e)
}, Ray.prototype.isIntersectionPlane = function(e) {
	return console.warn("THREE.Ray: .isIntersectionPlane() has been renamed to .intersectsPlane()."), this.intersectsPlane(e)
}, Ray.prototype.isIntersectionSphere = function(e) {
	return console.warn("THREE.Ray: .isIntersectionSphere() has been renamed to .intersectsSphere()."), this.intersectsSphere(e)
}, Triangle.prototype.area = function() {
	return console.warn("THREE.Triangle: .area() has been renamed to .getArea()."), this.getArea()
}, Triangle.prototype.barycoordFromPoint = function(e, t) {
	return console.warn("THREE.Triangle: .barycoordFromPoint() has been renamed to .getBarycoord()."), this.getBarycoord(e, t)
}, Triangle.prototype.midpoint = function(e) {
	return console.warn("THREE.Triangle: .midpoint() has been renamed to .getMidpoint()."), this.getMidpoint(e)
}, Triangle.prototypenormal = function(e) {
	return console.warn("THREE.Triangle: .normal() has been renamed to .getNormal()."), this.getNormal(e)
}, Triangle.prototype.plane = function(e) {
	return console.warn("THREE.Triangle: .plane() has been renamed to .getPlane()."), this.getPlane(e)
}, Triangle.barycoordFromPoint = function(e, t, i, n, r) {
	return console.warn("THREE.Triangle: .barycoordFromPoint() has been renamed to .getBarycoord()."), Triangle.getBarycoord(e, t, i, n, r)
}, Triangle.normal = function(e, t, i, n) {
	return console.warn("THREE.Triangle: .normal() has been renamed to .getNormal()."), Triangle.getNormal(e, t, i, n)
}, Vector2.prototype.fromAttribute = function(e, t, i) {
	return console.warn("THREE.Vector2: .fromAttribute() has been renamed to .fromBufferAttribute()."), this.fromBufferAttribute(e, t, i)
}, Vector2.prototype.distanceToManhattan = function(e) {
	return console.warn("THREE.Vector2: .distanceToManhattan() has been renamed to .manhattanDistanceTo()."), this.manhattanDistanceTo(e)
}, Vector2.prototype.lengthManhattan = function() {
	return console.warn("THREE.Vector2: .lengthManhattan() has been renamed to .manhattanLength()."), this.manhattanLength()
}, Vector3.prototype.setEulerFromRotationMatrix = function() {
	console.error("THREE.Vector3: .setEulerFromRotationMatrix() has been removed. Use Euler.setFromRotationMatrix() instead.")
}, Vector3.prototype.setEulerFromQuaternion = function() {
	console.error("THREE.Vector3: .setEulerFromQuaternion() has been removed. Use Euler.setFromQuaternion() instead.")
}, Vector3.prototype.getPositionFromMatrix = function(e) {
	return console.warn("THREE.Vector3: .getPositionFromMatrix() has been renamed to .setFromMatrixPosition()."), this.setFromMatrixPosition(e)
}, Vector3.prototype.getScaleFromMatrix = function(e) {
	return console.warn("THREE.Vector3: .getScaleFromMatrix() has been renamed to .setFromMatrixScale()."), this.setFromMatrixScale(e)
}, Vector3.prototype.getColumnFromMatrix = function(e, t) {
	return console.warn("THREE.Vector3: .getColumnFromMatrix() has been renamed to .setFromMatrixColumn()."), this.setFromMatrixColumn(t, e)
}, Vector3.prototype.applyProjection = function(e) {
	return console.warn("THREE.Vector3: .applyProjection() has been removed. Use .applyMatrix4( m ) instead."), this.applyMatrix4(e)
}, Vector3.prototype.fromAttribute = function(e, t, i) {
	return console.warn("THREE.Vector3: .fromAttribute() has been renamed to .fromBufferAttribute()."), this.fromBufferAttribute(e, t, i)
}, Vector3.prototype.distanceToManhattan = function(e) {
	return console.warn("THREE.Vector3: .distanceToManhattan() has been renamed to .manhattanDistanceTo()."), this.manhattanDistanceTo(e)
}, Vector3.prototype.lengthManhattan = function() {
	return console.warn("THREE.Vector3: .lengthManhattan() has been renamed to .manhattanLength()."), this.manhattanLength()
}, Vector4.prototype.fromAttribute = function(e, t, i) {
	return console.warn("THREE.Vector4: .fromAttribute() has been renamed to .fromBufferAttribute()."), this.fromBufferAttribute(e, t, i)
}, Vector4.prototype.lengthManhattan = function() {
	return console.warn("THREE.Vector4: .lengthManhattan() has been renamed to .manhattanLength()."), this.manhattanLength()
}, Object3D.prototype.getChildByName = function(e) {
	return console.warn("THREE.Object3D: .getChildByName() has been renamed to .getObjectByName()."), this.getObjectByName(e)
}, Object3D.prototype.renderDepth = function() {
	console.warn("THREE.Object3D: .renderDepth has been removed. Use .renderOrder, instead.")
}, Object3D.prototype.translate = function(e, t) {
	return console.warn("THREE.Object3D: .translate() has been removed. Use .translateOnAxis( axis, distance ) instead."), this.translateOnAxis(t, e)
}, Object3D.prototype.getWorldRotation = function() {
	console.error("THREE.Object3D: .getWorldRotation() has been removed. Use THREE.Object3D.getWorldQuaternion( target ) instead.")
}, Object3D.prototype.applyMatrix = function(e) {
	return console.warn("THREE.Object3D: .applyMatrix() has been renamed to .applyMatrix4()."), this.applyMatrix4(e)
}, Object.defineProperties(Object3D.prototype, {
	eulerOrder: {
		get: function() {
			return console.warn("THREE.Object3D: .eulerOrder is now .rotation.order."), this.rotation.order
		},
		set: function(e) {
			console.warn("THREE.Object3D: .eulerOrder is now .rotation.order."), this.rotation.order = e
		}
	},
	useQuaternion: {
		get: function() {
			console.warn("THREE.Object3D: .useQuaternion has been removed. The library now uses quaternions by default.")
		},
		set: function() {
			console.warn("THREE.Object3D: .useQuaternion has been removed. The library now uses quaternions by default.")
		}
	}
}), Mesh.prototype.setDrawMode = function() {
	console.error("THREE.Mesh: .setDrawMode() has been removed. The renderer now always assumes THREE.TrianglesDrawMode. Transform your geometry via BufferGeometryUtils.toTrianglesDrawMode() if necessary.")
}, Object.defineProperties(Mesh.prototype, {
	drawMode: {
		get: function() {
			return console.error("THREE.Mesh: .drawMode has been removed. The renderer now always assumes THREE.TrianglesDrawMode."), 0
		},
		set: function() {
			console.error("THREE.Mesh: .drawMode has been removed. The renderer now always assumes THREE.TrianglesDrawMode. Transform your geometry via BufferGeometryUtils.toTrianglesDrawMode() if necessary.")
		}
	}
}), Object.defineProperties(LOD.prototype, {
	objects: {
		get: function() {
			return console.warn("THREE.LOD: .objects has been renamed to .levels."), this.levels
		}
	}
}), Object.defineProperty(Skeleton.prototype, "useVertexTexture", {
	get: function() {
		console.warn("THREE.Skeleton: useVertexTexture has been removed.")
	},
	set: function() {
		console.warn("THREE.Skeleton: useVertexTexture has been removed.")
	}
}), SkinnedMesh.prototype.initBones = function() {
	console.error("THREE.SkinnedMesh: initBones() has been removed.")
}, Object.defineProperty(Curve.prototype, "__arcLengthDivisions", {
	get: function() {
		return console.warn("THREE.Curve: .__arcLengthDivisions is now .arcLengthDivisions."), this.arcLengthDivisions
	},
	set: function(e) {
		console.warn("THREE.Curve: .__arcLengthDivisions is now .arcLengthDivisions."), this.arcLengthDivisions = e
	}
}), PerspectiveCamera.prototype.setLens = function(e, t) {
	console.warn("THREE.PerspectiveCamera.setLens is deprecated. Use .setFocalLength and .filmGauge for a photographic setup."), void 0 !== t && (this.filmGauge = t), this.setFocalLength(e)
}, Object.defineProperties(Light.prototype, {
	onlyShadow: {
		set: function() {
			console.warn("THREE.Light: .onlyShadow has been removed.")
		}
	},
	shadowCameraFov: {
		set: function(e) {
			console.warn("THREE.Light: .shadowCameraFov is now .shadow.camera.fov."), this.shadow.camera.fov = e
		}
	},
	shadowCameraLeft: {
		set: function(e) {
			console.warn("THREE.Light: .shadowCameraLeft is now .shadow.camera.left."), this.shadow.camera.left = e
		}
	},
	shadowCameraRight: {
		set: function(e) {
			console.warn("THREE.Light: .shadowCameraRight is now .shadow.camera.right."), this.shadow.camera.right = e
		}
	},
	shadowCameraTop: {
		set: function(e) {
			console.warn("THREE.Light: .shadowCameraTop is now .shadow.camera.top."), this.shadow.camera.top = e
		}
	},
	shadowCameraBottom: {
		set: function(e) {
			console.warn("THREE.Light: .shadowCameraBottom is now .shadow.camera.bottom."), this.shadow.camera.bottom = e
		}
	},
	shadowCameraNear: {
		set: function(e) {
			console.warn("THREE.Light: .shadowCameraNear is now .shadow.camera.near."), this.shadow.camera.near = e
		}
	},
	shadowCameraFar: {
		set: function(e) {
			console.warn("THREE.Light: .shadowCameraFar is now .shadow.camera.far."), this.shadow.camera.far = e
		}
	},
	shadowCameraVisible: {
		set: function() {
			console.warn("THREE.Light: .shadowCameraVisible has been removed. Use new THREE.CameraHelper( light.shadow.camera ) instead.")
		}
	},
	shadowBias: {
		set: function(e) {
			console.warn("THREE.Light: .shadowBias is now .shadow.bias."), this.shadow.bias = e
		}
	},
	shadowDarkness: {
		set: function() {
			console.warn("THREE.Light: .shadowDarkness has been removed.")
		}
	},
	shadowMapWidth: {
		set: function(e) {
			console.warn("THREE.Light: .shadowMapWidth is now .shadow.mapSize.width."), this.shadow.mapSize.width = e
		}
	},
	shadowMapHeight: {
		set: function(e) {
			console.warn("THREE.Light: .shadowMapHeight is now .shadow.mapSize.height."), this.shadow.mapSize.height = e
		}
	}
}), Object.defineProperties(BufferAttribute.prototype, {
	length: {
		get: function() {
			return console.warn("THREE.BufferAttribute: .length has been deprecated. Use .count instead."), this.array.length
		}
	},
	dynamic: {
		get: function() {
			return console.warn("THREE.BufferAttribute: .dynamic has been deprecated. Use .usage instead."), 35048 === this.usage
		},
		set: function() {
			console.warn("THREE.BufferAttribute: .dynamic has been deprecated. Use .usage instead."), this.setUsage(35048)
		}
	}
}), BufferAttribute.prototype.setDynamic = function(e) {
	return console.warn("THREE.BufferAttribute: .setDynamic() has been deprecated. Use .setUsage() instead."), this.setUsage(!0 === e ? 35048 : 35044), this
}, BufferAttribute.prototype.copyIndicesArray = function() {
	console.error("THREE.BufferAttribute: .copyIndicesArray() has been removed.")
}, BufferAttribute.prototype.setArray = function() {
	console.error("THREE.BufferAttribute: .setArray has been removed. Use BufferGeometry .setAttribute to replace/resize attribute buffers")
}, BufferGeometry.prototype.addIndex = function(e) {
	console.warn("THREE.BufferGeometry: .addIndex() has been renamed to .setIndex()."), this.setIndex(e)
}, BufferGeometry.prototype.addAttribute = function(e, t) {
	return console.warn("THREE.BufferGeometry: .addAttribute() has been renamed to .setAttribute()."), t && t.isBufferAttribute || t && t.isInterleavedBufferAttribute ? "index" === e ? (console.warn("THREE.BufferGeometry.addAttribute: Use .setIndex() for index attribute."), this.setIndex(t), this) : this.setAttribute(e, t) : (console.warn("THREE.BufferGeometry: .addAttribute() now expects ( name, attribute )."), this.setAttribute(e, new BufferAttribute(arguments[1], arguments[2])))
}, BufferGeometry.prototype.addDrawCall = function(e, t, i) {
	void 0 !== i && console.warn("THREE.BufferGeometry: .addDrawCall() no longer supports indexOffset."), console.warn("THREE.BufferGeometry: .addDrawCall() is now .addGroup()."), this.addGroup(e, t)
}, BufferGeometry.prototype.clearDrawCalls = function() {
	console.warn("THREE.BufferGeometry: .clearDrawCalls() is now .clearGroups()."), this.clearGroups()
}, BufferGeometry.prototype.computeOffsets = function() {
	console.warn("THREE.BufferGeometry: .computeOffsets() has been removed.")
}, BufferGeometry.prototype.removeAttribute = function(e) {
	return console.warn("THREE.BufferGeometry: .removeAttribute() has been renamed to .deleteAttribute()."), this.deleteAttribute(e)
}, BufferGeometry.prototype.applyMatrix = function(e) {
	return console.warn("THREE.BufferGeometry: .applyMatrix() has been renamed to .applyMatrix4()."), this.applyMatrix4(e)
}, Object.defineProperties(BufferGeometry.prototype, {
	drawcalls: {
		get: function() {
			return console.error("THREE.BufferGeometry: .drawcalls has been renamed to .groups."), this.groups
		}
	},
	offsets: {
		get: function() {
			return console.warn("THREE.BufferGeometry: .offsets has been renamed to .groups."), this.groups
		}
	}
}), Object.defineProperties(InstancedBufferGeometry.prototype, {
	maxInstancedCount: {
		get: function() {
			return console.warn("THREE.InstancedBufferGeometry: .maxInstancedCount has been renamed to .instanceCount."), this.instanceCount
		},
		set: function(e) {
			console.warn("THREE.InstancedBufferGeometry: .maxInstancedCount has been renamed to .instanceCount."), this.instanceCount = e
		}
	}
}), Object.defineProperties(Raycaster.prototype, {
	linePrecision: {
		get: function() {
			return console.warn("THREE.Raycaster: .linePrecision has been deprecated. Use .params.Line.threshold instead."), this.params.Line.threshold
		},
		set: function(e) {
			console.warn("THREE.Raycaster: .linePrecision has been deprecated. Use .params.Line.threshold instead."), this.params.Line.threshold = e
		}
	}
}), Object.defineProperties(InterleavedBuffer.prototype, {
	dynamic: {
		get: function() {
			return console.warn("THREE.InterleavedBuffer: .length has been deprecated. Use .usage instead."), 35048 === this.usage
		},
		set: function(e) {
			console.warn("THREE.InterleavedBuffer: .length has been deprecated. Use .usage instead."), this.setUsage(e)
		}
	}
}), InterleavedBuffer.prototype.setDynamic = function(e) {
	return console.warn("THREE.InterleavedBuffer: .setDynamic() has been deprecated. Use .setUsage() instead."), this.setUsage(!0 === e ? 35048 : 35044), this
}, InterleavedBuffer.prototype.setArray = function() {
	console.error("THREE.InterleavedBuffer: .setArray has been removed. Use BufferGeometry .setAttribute to replace/resize attribute buffers")
}, Scene.prototype.dispose = function() {
	console.error("THREE.Scene: .dispose() has been removed.")
}, Object.defineProperties(Uniform.prototype, {
	dynamic: {
		set: function() {
			console.warn("THREE.Uniform: .dynamic has been removed. Use object.onBeforeRender() instead.")
		}
	},
	onUpdate: {
		value: function() {
			return console.warn("THREE.Uniform: .onUpdate() has been removed. Use object.onBeforeRender() instead."), this
		}
	}
}), Object.defineProperties(Material$1.prototype, {
	wrapAround: {
		get: function() {
			console.warn("THREE.Material: .wrapAround has been removed.")
		},
		set: function() {
			console.warn("THREE.Material: .wrapAround has been removed.")
		}
	},
	overdraw: {
		get: function() {
			console.warn("THREE.Material: .overdraw has been removed.")
		},
		set: function() {
			console.warn("THREE.Material: .overdraw has been removed.")
		}
	},
	wrapRGB: {
		get: function() {
			return console.warn("THREE.Material: .wrapRGB has been removed."), new Color
		}
	},
	shading: {
		get: function() {
			console.error("THREE." + this.type + ": .shading has been removed. Use the boolean .flatShading instead.")
		},
		set: function(e) {
			console.warn("THREE." + this.type + ": .shading has been removed. Use the boolean .flatShading instead."), this.flatShading = 1 === e
		}
	},
	stencilMask: {
		get: function() {
			return console.warn("THREE." + this.type + ": .stencilMask has been removed. Use .stencilFuncMask instead."), this.stencilFuncMask
		},
		set: function(e) {
			console.warn("THREE." + this.type + ": .stencilMask has been removed. Use .stencilFuncMask instead."), this.stencilFuncMask = e
		}
	}
}), Object.defineProperties(MeshPhongMaterial.prototype, {
	metal: {
		get: function() {
			return console.warn("THREE.MeshPhongMaterial: .metal has been removed. Use THREE.MeshStandardMaterial instead."), !1
		},
		set: function() {
			console.warn("THREE.MeshPhongMaterial: .metal has been removed. Use THREE.MeshStandardMaterial instead")
		}
	}
}), Object.defineProperties(MeshPhysicalMaterial.prototype, {
	transparency: {
		get: function() {
			return console.warn("THREE.MeshPhysicalMaterial: .transparency has been renamed to .transmission."), this.transmission
		},
		set: function(e) {
			console.warn("THREE.MeshPhysicalMaterial: .transparency has been renamed to .transmission."), this.transmission = e
		}
	}
}), Object.defineProperties(ShaderMaterial.prototype, {
	derivatives: {
		get: function() {
			return console.warn("THREE.ShaderMaterial: .derivatives has been moved to .extensions.derivatives."), this.extensions.derivatives
		},
		set: function(e) {
			console.warn("THREE. ShaderMaterial: .derivatives has been moved to .extensions.derivatives."), this.extensions.derivatives = e
		}
	}
}), WebGLRenderer.prototype.clearTarget = function(e, t, i, n) {
	console.warn("THREE.WebGLRenderer: .clearTarget() has been deprecated. Use .setRenderTarget() and .clear() instead."), this.setRenderTarget(e), this.clear(t, i, n)
}, WebGLRenderer.prototype.animate = function(e) {
	console.warn("THREE.WebGLRenderer: .animate() is now .setAnimationLoop()."), this.setAnimationLoop(e)
}, WebGLRenderer.prototype.getCurrentRenderTarget = function() {
	return console.warn("THREE.WebGLRenderer: .getCurrentRenderTarget() is now .getRenderTarget()."), this.getRenderTarget()
}, WebGLRenderer.prototype.getMaxAnisotropy = function() {
	return console.warn("THREE.WebGLRenderer: .getMaxAnisotropy() is now .capabilities.getMaxAnisotropy()."), this.capabilities.getMaxAnisotropy()
}, WebGLRenderer.prototype.getPrecision = function() {
	return console.warn("THREE.WebGLRenderer: .getPrecision() is now .capabilities.precision."), this.capabilities.precision
}, WebGLRenderer.prototype.resetGLState = function() {
	return console.warn("THREE.WebGLRenderer: .resetGLState() is now .state.reset()."), this.state.reset()
}, WebGLRenderer.prototype.supportsFloatTextures = function() {
	return console.warn("THREE.WebGLRenderer: .supportsFloatTextures() is now .extensions.get( 'OES_texture_float' )."), this.extensions.get("OES_texture_float")
}, WebGLRenderer.prototype.supportsHalfFloatTextures = function() {
	return console.warn("THREE.WebGLRenderer: .supportsHalfFloatTextures() is now .extensions.get( 'OES_texture_half_float' )."), this.extensions.get("OES_texture_half_float")
}, WebGLRenderer.prototype.supportsStandardDerivatives = function() {
	return console.warn("THREE.WebGLRenderer: .supportsStandardDerivatives() is now .extensions.get( 'OES_standard_derivatives' )."), this.extensions.get("OES_standard_derivatives")
}, WebGLRenderer.prototype.supportsCompressedTextureS3TC = function() {
	return console.warn("THREE.WebGLRenderer: .supportsCompressedTextureS3TC() is now .extensions.get( 'WEBGL_compressed_texture_s3tc' )."), this.extensions.get("WEBGL_compressed_texture_s3tc")
}, WebGLRenderer.prototype.supportsCompressedTexturePVRTC = function() {
	return console.warn("THREE.WebGLRenderer: .supportsCompressedTexturePVRTC() is now .extensions.get( 'WEBGL_compressed_texture_pvrtc' )."), this.extensions.get("WEBGL_compressed_texture_pvrtc")
}, WebGLRenderer.prototype.supportsBlendMinMax = function() {
	return console.warn("THREE.WebGLRenderer: .supportsBlendMinMax() is now .extensions.get( 'EXT_blend_minmax' )."), this.extensions.get("EXT_blend_minmax")
}, WebGLRenderer.prototype.supportsVertexTextures = function() {
	return console.warn("THREE.WebGLRenderer: .supportsVertexTextures() is now .capabilities.vertexTextures."), this.capabilities.vertexTextures
}, WebGLRenderer.prototype.supportsInstancedArrays = function() {
	return console.warn("THREE.WebGLRenderer: .supportsInstancedArrays() is now .extensions.get( 'ANGLE_instanced_arrays' )."), this.extensions.get("ANGLE_instanced_arrays")
}, WebGLRenderer.prototype.enableScissorTest = function(e) {
	console.warn("THREE.WebGLRenderer: .enableScissorTest() is now .setScissorTest()."), this.setScissorTest(e)
}, WebGLRenderer.prototype.initMaterial = function() {
	console.warn("THREE.WebGLRenderer: .initMaterial() has been removed.")
}, WebGLRenderer.prototype.addPrePlugin = function() {
	console.warn("THREE.WebGLRenderer: .addPrePlugin() has been removed.")
}, WebGLRenderer.prototype.addPostPlugin = function() {
	console.warn("THREE.WebGLRenderer: .addPostPlugin() has been removed.")
}, WebGLRenderer.prototype.updateShadowMap = function() {
	console.warn("THREE.WebGLRenderer: .updateShadowMap() has been removed.")
}, WebGLRenderer.prototype.setFaceCulling = function() {
	console.warn("THREE.WebGLRenderer: .setFaceCulling() has been removed.")
}, WebGLRenderer.prototype.allocTextureUnit = function() {
	console.warn("THREE.WebGLRenderer: .allocTextureUnit() has been removed.")
}, WebGLRenderer.prototype.setTexture = function() {
	console.warn("THREE.WebGLRenderer: .setTexture() has been removed.")
}, WebGLRenderer.prototype.setTexture2D = function() {
	console.warn("THREE.WebGLRenderer: .setTexture2D() has been removed.")
}, WebGLRenderer.prototype.setTextureCube = function() {
	console.warn("THREE.WebGLRenderer: .setTextureCube() has been removed.")
}, WebGLRenderer.prototype.getActiveMipMapLevel = function() {
	return console.warn("THREE.WebGLRenderer: .getActiveMipMapLevel() is now .getActiveMipmapLevel()."), this.getActiveMipmapLevel()
}, Object.defineProperties(WebGLRenderer.prototype, {
	shadowMapEnabled: {
		get: function() {
			return this.shadowMap.enabled
		},
		set: function(e) {
			console.warn("THREE.WebGLRenderer: .shadowMapEnabled is now .shadowMap.enabled."), this.shadowMap.enabled = e
		}
	},
	shadowMapType: {
		get: function() {
			return this.shadowMap.type
		},
		set: function(e) {
			console.warn("THREE.WebGLRenderer: .shadowMapType is now .shadowMap.type."), this.shadowMap.type = e
		}
	},
	shadowMapCullFace: {
		get: function() {
			console.warn("THREE.WebGLRenderer: .shadowMapCullFace has been removed. Set Material.shadowSide instead.")
		},
		set: function() {
			console.warn("THREE.WebGLRenderer: .shadowMapCullFace has been removed. Set Material.shadowSide instead.")
		}
	},
	context: {
		get: function() {
			return console.warn("THREE.WebGLRenderer: .context has been removed. Use .getContext() instead."), this.getContext()
		}
	},
	vr: {
		get: function() {
			return console.warn("THREE.WebGLRenderer: .vr has been renamed to .xr"), this.xr
		}
	},
	gammaInput: {
		get: function() {
			return console.warn("THREE.WebGLRenderer: .gammaInput has been removed. Set the encoding for textures via Texture.encoding instead."), !1
		},
		set: function() {
			console.warn("THREE.WebGLRenderer: .gammaInput has been removed. Set the encoding for textures via Texture.encoding instead.")
		}
	},
	gammaOutput: {
		get: function() {
			return console.warn("THREE.WebGLRenderer: .gammaOutput has been removed. Set WebGLRenderer.outputEncoding instead."), !1
		},
		set: function(e) {
			console.warn("THREE.WebGLRenderer: .gammaOutput has been removed. Set WebGLRenderer.outputEncoding instead."), this.outputEncoding = !0 === e ? 3001 : 3e3
		}
	},
	toneMappingWhitePoint: {
		get: function() {
			return console.warn("THREE.WebGLRenderer: .toneMappingWhitePoint has been removed."), 1
		},
		set: function() {
			console.warn("THREE.WebGLRenderer: .toneMappingWhitePoint has been removed.")
		}
	}
}), Object.defineProperties(WebGLShadowMap.prototype, {
	cullFace: {
		get: function() {
			console.warn("THREE.WebGLRenderer: .shadowMap.cullFace has been removed. Set Material.shadowSide instead.")
		},
		set: function() {
			console.warn("THREE.WebGLRenderer: .shadowMap.cullFace has been removed. Set Material.shadowSide instead.")
		}
	},
	renderReverseSided: {
		get: function() {
			console.warn("THREE.WebGLRenderer: .shadowMap.renderReverseSided has been removed. Set Material.shadowSide instead.")
		},
		set: function() {
			console.warn("THREE.WebGLRenderer: .shadowMap.renderReverseSided has been removed. Set Material.shadowSide instead.")
		}
	},
	renderSingleSided: {
		get: function() {
			console.warn("THREE.WebGLRenderer: .shadowMap.renderSingleSided has been removed. Set Material.shadowSide instead.")
		},
		set: function() {
			console.warn("THREE.WebGLRenderer: .shadowMap.renderSingleSided has been removed. Set Material.shadowSide instead.")
		}
	}
}), Object.defineProperties(WebGLRenderTarget.prototype, {
	wrapS: {
		get: function() {
			return console.warn("THREE.WebGLRenderTarget: .wrapS is now .texture.wrapS."), this.texture.wrapS
		},
		set: function(e) {
			console.warn("THREE.WebGLRenderTarget: .wrapS is now .texture.wrapS."), this.texture.wrapS = e
		}
	},
	wrapT: {
		get: function() {
			return console.warn("THREE.WebGLRenderTarget: .wrapT is now .texture.wrapT."), this.texture.wrapT
		},
		set: function(e) {
			console.warn("THREE.WebGLRenderTarget: .wrapT is now .texture.wrapT."), this.texture.wrapT = e
		}
	},
	magFilter: {
		get: function() {
			return console.warn("THREE.WebGLRenderTarget: .magFilter is now .texture.magFilter."), this.texture.magFilter
		},
		set: function(e) {
			console.warn("THREE.WebGLRenderTarget: .magFilter is now .texture.magFilter."), this.texture.magFilter = e
		}
	},
	minFilter: {
		get: function() {
			return console.warn("THREE.WebGLRenderTarget: .minFilter is now .texture.minFilter."), this.texture.minFilter
		},
		set: function(e) {
			console.warn("THREE.WebGLRenderTarget: .minFilter is now .texture.minFilter."), this.texture.minFilter = e
		}
	},
	anisotropy: {
		get: function() {
			return console.warn("THREE.WebGLRenderTarget: .anisotropy is now .texture.anisotropy."), this.texture.anisotropy
		},
		set: function(e) {
			console.warn("THREE.WebGLRenderTarget: .anisotropy is now .texture.anisotropy."), this.texture.anisotropy = e
		}
	},
	offset: {
		get: function() {
			return console.warn("THREE.WebGLRenderTarget: .offset is now .texture.offset."), this.texture.offset
		},
		set: function(e) {
			console.warn("THREE.WebGLRenderTarget: .offset is now .texture.offset."), this.texture.offset = e
		}
	},
	repeat: {
		get: function() {
			return console.warn("THREE.WebGLRenderTarget: .repeat is now .texture.repeat."), this.texture.repeat
		},
		set: function(e) {
			console.warn("THREE.WebGLRenderTarget: .repeat is now .texture.repeat."), this.texture.repeat = e
		}
	},
	format: {
		get: function() {
			return console.warn("THREE.WebGLRenderTarget: .format is now .texture.format."), this.texture.format
		},
		set: function(e) {
			console.warn("THREE.WebGLRenderTarget: .format is now .texture.format."), this.texture.format = e
		}
	},
	type: {
		get: function() {
			return console.warn("THREE.WebGLRenderTarget: .type is now .texture.type."), this.texture.type
		},
		set: function(e) {
			console.warn("THREE.WebGLRenderTarget: .type is now .texture.type."), this.texture.type = e
		}
	},
	generateMipmaps: {
		get: function() {
			return console.warn("THREE.WebGLRenderTarget: .generateMipmaps is now .texture.generateMipmaps."), this.texture.generateMipmaps
		},
		set: function(e) {
			console.warn("THREE.WebGLRenderTarget: .generateMipmaps is now .texture.generateMipmaps."), this.texture.generateMipmaps = e
		}
	}
}), Object.defineProperties(Audio.prototype, {
	load: {
		value: function(e) {
			console.warn("THREE.Audio: .load has been deprecated. Use THREE.AudioLoader instead.");
			const t = this;
			return (new AudioLoader).load(e, (function(e) {
				t.setBuffer(e)
			})), this
		}
	},
	startTime: {
		set: function() {
			console.warn("THREE.Audio: .startTime is now .play( delay ).")
		}
	}
}), CubeCamera.prototype.updateCubeMap = function(e, t) {
	return console.warn("THREE.CubeCamera: .updateCubeMap() is now .update()."), this.update(e, t)
}, CubeCamera.prototype.clear = function(e, t, i, n) {
	return console.warn("THREE.CubeCamera: .clear() is now .renderTarget.clear()."), this.renderTarget.clear(e, t, i, n)
}, ImageUtils.crossOrigin = void 0, ImageUtils.loadTexture = function(e, t, i, n) {
	console.warn("THREE.ImageUtils.loadTexture has been deprecated. Use THREE.TextureLoader() instead.");
	const r = new TextureLoader;
	r.setCrossOrigin(this.crossOrigin);
	const a = r.load(e, i, void 0, n);
	return t && (a.mapping = t), a
}, ImageUtils.loadTextureCube = function(e, t, i, n) {
	console.warn("THREE.ImageUtils.loadTextureCube has been deprecated. Use THREE.CubeTextureLoader() instead.");
	const r = new CubeTextureLoader;
	r.setCrossOrigin(this.crossOrigin);
	const a = r.load(e, i, void 0, n);
	return t && (a.mapping = t), a
}, ImageUtils.loadCompressedTexture = function() {
	console.error("THREE.ImageUtils.loadCompressedTexture has been removed. Use THREE.DDSLoader instead.")
}, ImageUtils.loadCompressedTextureCube = function() {
	console.error("THREE.ImageUtils.loadCompressedTextureCube has been removed. Use THREE.DDSLoader instead.")
}, "undefined" != typeof __THREE_DEVTOOLS__ && __THREE_DEVTOOLS__.dispatchEvent(new CustomEvent("register", {
	detail: {
		revision: "126"
	}
})), "undefined" != typeof window && (window.__THREE__ ? console.warn("WARNING: Multiple instances of Three.js being imported.") : window.__THREE__ = "126");
var DRACOLoader = function(e) {
	Loader.call(this, e), this.decoderPath = "", this.decoderConfig = {}, this.decoderBinary = null, this.decoderPending = null, this.workerLimit = 4, this.workerPool = [], this.workerNextTaskID = 1, this.workerSourceURL = "", this.defaultAttributeIDs = {
		position: "POSITION",
		normal: "NORMAL",
		color: "COLOR",
		uv: "TEX_COORD"
	}, this.defaultAttributeTypes = {
		position: "Float32Array",
		normal: "Float32Array",
		color: "Float32Array",
		uv: "Float32Array"
	}
};
DRACOLoader.prototype = Object.assign(Object.create(Loader.prototype), {
	constructor: DRACOLoader,
	setDecoderPath: function(e) {
		return this.decoderPath = e, this
	},
	setDecoderConfig: function(e) {
		return this.decoderConfig = e, this
	},
	setWorkerLimit: function(e) {
		return this.workerLimit = e, this
	},
	setVerbosity: function() {
		console.warn("THREE.DRACOLoader: The .setVerbosity() method has been removed.")
	},
	setDrawMode: function() {
		console.warn("THREE.DRACOLoader: The .setDrawMode() method has been removed.")
	},
	setSkipDequantization: function() {
		console.warn("THREE.DRACOLoader: The .setSkipDequantization() method has been removed.")
	},
	load: function(e, t, i, n) {
		var r = new FileLoader(this.manager);
		r.setPath(this.path), r.setResponseType("arraybuffer"), r.setRequestHeader(this.requestHeader), r.setWithCredentials(this.withCredentials), r.load(e, e => {
			var i = {
				attributeIDs: this.defaultAttributeIDs,
				attributeTypes: this.defaultAttributeTypes,
				useUniqueIDs: !1
			};
			this.decodeGeometry(e, i).then(t).catch(n)
		}, i, n)
	},
	decodeDracoFile: function(e, t, i, n) {
		var r = {
			attributeIDs: i || this.defaultAttributeIDs,
			attributeTypes: n || this.defaultAttributeTypes,
			useUniqueIDs: !!i
		};
		this.decodeGeometry(e, r).then(t)
	},
	decodeGeometry: function(e, t) {
		for (var i in t.attributeTypes) {
			var n = t.attributeTypes[i];
			void 0 !== n.BYTES_PER_ELEMENT && (t.attributeTypes[i] = n.name)
		}
		var r, a = JSON.stringify(t);
		if (DRACOLoader.taskCache.has(e)) {
			var s = DRACOLoader.taskCache.get(e);
			if (s.key === a) return s.promise;
			if (0 === e.byteLength) throw new Error("THREE.DRACOLoader: Unable to re-decode a buffer with different settings. Buffer has already been transferred.")
		}
		var o = this.workerNextTaskID++,
			l = e.byteLength,
			c = this._getWorker(o, l).then(i => (r = i, new Promise((i, n) => {
				r._callbacks[o] = {
					resolve: i,
					reject: n
				}, r.postMessage({
					type: "decode",
					id: o,
					taskConfig: t,
					buffer: e
				}, [e])
			}))).then(e => this._createGeometry(e.geometry));
		return c.catch(() => !0).then(() => {
			r && o && this._releaseTask(r, o)
		}), DRACOLoader.taskCache.set(e, {
			key: a,
			promise: c
		}), c
	},
	_createGeometry: function(e) {
		var t = new BufferGeometry;
		e.index && t.setIndex(new BufferAttribute(e.index.array, 1));
		for (var i = 0; i < e.attributes.length; i++) {
			var n = e.attributes[i],
				r = n.name,
				a = n.array,
				s = n.itemSize;
			t.setAttribute(r, new BufferAttribute(a, s))
		}
		return t
	},
	_loadLibrary: function(e, t) {
		var i = new FileLoader(this.manager);
		return i.setPath(this.decoderPath), i.setResponseType(t), i.setWithCredentials(this.withCredentials), new Promise((t, n) => {
			i.load(e, t, void 0, n)
		})
	},
	preload: function() {
		return this._initDecoder(), this
	},
	_initDecoder: function() {
		if (this.decoderPending) return this.decoderPending;
		var e = "object" != typeof WebAssembly || "js" === this.decoderConfig.type,
			t = [];
		return e ? t.push(this._loadLibrary("draco_decoder.js", "text")) : (t.push(this._loadLibrary("draco_wasm_wrapper.js", "text")), t.push(this._loadLibrary("draco_decoder.wasm", "arraybuffer"))), this.decoderPending = Promise.all(t).then(t => {
			var i = t[0];
			e || (this.decoderConfig.wasmBinary = t[1]);
			var n = DRACOLoader.DRACOWorker.toString(),
				r = ["/* draco decoder */", i, "", "/* worker */", n.substring(n.indexOf("{") + 1, n.lastIndexOf("}"))].join("\n");
			this.workerSourceURL = URL.createObjectURL(new Blob([r]))
		}), this.decoderPending
	},
	_getWorker: function(e, t) {
		return this._initDecoder().then(() => {
			var i;
			this.workerPool.length < this.workerLimit ? ((i = new Worker(this.workerSourceURL))._callbacks = {}, i._taskCosts = {}, i._taskLoad = 0, i.postMessage({
				type: "init",
				decoderConfig: this.decoderConfig
			}), i.onmessage = function(e) {
				var t = e.data;
				switch (t.type) {
					case "decode":
						i._callbacks[t.id].resolve(t);
						break;
					case "error":
						i._callbacks[t.id].reject(t);
						break;
					default:
						console.error('THREE.DRACOLoader: Unexpected message, "' + t.type + '"')
				}
			}, this.workerPool.push(i)) : this.workerPool.sort((function(e, t) {
				return e._taskLoad > t._taskLoad ? -1 : 1
			}));
			return (i = this.workerPool[this.workerPool.length - 1])._taskCosts[e] = t, i._taskLoad += t, i
		})
	},
	_releaseTask: function(e, t) {
		e._taskLoad -= e._taskCosts[t], delete e._callbacks[t], delete e._taskCosts[t]
	},
	debug: function() {
		console.log("Task load: ", this.workerPool.map(e => e._taskLoad))
	},
	dispose: function() {
		for (var e = 0; e < this.workerPool.length; ++e) this.workerPool[e].terminate();
		return this.workerPool.length = 0, this
	}
}), DRACOLoader.DRACOWorker = function() {
	var e, t;

	function i(e, t, i, n, r, a) {
		var s = a.num_components(),
			o = i.num_points() * s,
			l = o * r.BYTES_PER_ELEMENT,
			c = function(e, t) {
				switch (t) {
					case Float32Array:
						return e.DT_FLOAT32;
					case Int8Array:
						return e.DT_INT8;
					case Int16Array:
						return e.DT_INT16;
					case Int32Array:
						return e.DT_INT32;
					case Uint8Array:
						return e.DT_UINT8;
					case Uint16Array:
						return e.DT_UINT16;
					case Uint32Array:
						return e.DT_UINT32
				}
			}(e, r),
			h = e._malloc(l);
		t.GetAttributeDataArrayForAllPoints(i, a, c, l, h);
		var u = new r(e.HEAPF32.buffer, h, o).slice();
		return e._free(h), {
			name: n,
			array: u,
			itemSize: s
		}
	}
	onmessage = function(n) {
		var r = n.data;
		switch (r.type) {
			case "init":
				e = r.decoderConfig, t = new Promise((function(t) {
					e.onModuleLoaded = function(e) {
						t({
							draco: e
						})
					}, DracoDecoderModule(e)
				}));
				break;
			case "decode":
				var a = r.buffer,
					s = r.taskConfig;
				t.then(e => {
					var t = e.draco,
						n = new t.Decoder,
						o = new t.DecoderBuffer;
					o.Init(new Int8Array(a), a.byteLength);
					try {
						var l = function(e, t, n, r) {
								var a, s, o = r.attributeIDs,
									l = r.attributeTypes,
									c = t.GetEncodedGeometryType(n);
								if (c === e.TRIANGULAR_MESH) a = new e.Mesh, s = t.DecodeBufferToMesh(n, a);
								else {
									if (c !== e.POINT_CLOUD) throw new Error("THREE.DRACOLoader: Unexpected geometry type.");
									a = new e.PointCloud, s = t.DecodeBufferToPointCloud(n, a)
								}
								if (!s.ok() || 0 === a.ptr) throw new Error("THREE.DRACOLoader: Decoding failed: " + s.error_msg());
								var h = {
									index: null,
									attributes: []
								};
								for (var u in o) {
									var d, p, m = self[l[u]];
									if (r.useUniqueIDs) p = o[u], d = t.GetAttributeByUniqueId(a, p);
									else {
										if (-1 === (p = t.GetAttributeId(a, e[o[u]]))) continue;
										d = t.GetAttribute(a, p)
									}
									h.attributes.push(i(e, t, a, u, m, d))
								}
								c === e.TRIANGULAR_MESH && (h.index = function(e, t, i) {
									var n = 3 * i.num_faces(),
										r = 4 * n,
										a = e._malloc(r);
									t.GetTrianglesUInt32Array(i, r, a);
									var s = new Uint32Array(e.HEAPF32.buffer, a, n).slice();
									return e._free(a), {
										array: s,
										itemSize: 1
									}
								}(e, t, a));
								return e.destroy(a), h
							}(t, n, o, s),
							c = l.attributes.map(e => e.array.buffer);
						l.index && c.push(l.index.array.buffer), self.postMessage({
							type: "decode",
							id: r.id,
							geometry: l
						}, c)
					} catch (e) {
						console.error(e), self.postMessage({
							type: "error",
							id: r.id,
							error: e.message
						})
					} finally {
						t.destroy(o), t.destroy(n)
					}
				})
		}
	}
}, DRACOLoader.taskCache = new WeakMap, DRACOLoader.setDecoderPath = function() {
	console.warn("THREE.DRACOLoader: The .setDecoderPath() method has been removed. Use instance methods.")
}, DRACOLoader.setDecoderConfig = function() {
	console.warn("THREE.DRACOLoader: The .setDecoderConfig() method has been removed. Use instance methods.")
}, DRACOLoader.releaseDecoderModule = function() {
	console.warn("THREE.DRACOLoader: The .releaseDecoderModule() method has been removed. Use instance methods.")
}, DRACOLoader.getDecoderModule = function() {
	console.warn("THREE.DRACOLoader: The .getDecoderModule() method has been removed. Use instance methods.")
};
var GLTFLoader = function() {
		function e(e) {
			Loader.call(this, e), this.dracoLoader = null, this.ktx2Loader = null, this.meshoptDecoder = null, this.pluginCallbacks = [], this.register((function(e) {
				return new a(e)
			})), this.register((function(e) {
				return new o(e)
			})), this.register((function(e) {
				return new l(e)
			})), this.register((function(e) {
				return new s(e)
			})), this.register((function(e) {
				return new n(e)
			})), this.register((function(e) {
				return new c(e)
			}))
		}

		function t() {
			var e = {};
			return {
				get: function(t) {
					return e[t]
				},
				add: function(t, i) {
					e[t] = i
				},
				remove: function(t) {
					delete e[t]
				},
				removeAll: function() {
					e = {}
				}
			}
		}
		e.prototype = Object.assign(Object.create(Loader.prototype), {
			constructor: e,
			load: function(e, t, i, n) {
				var r, a = this;
				r = "" !== this.resourcePath ? this.resourcePath : "" !== this.path ? this.path : LoaderUtils.extractUrlBase(e), this.manager.itemStart(e);
				var s = function(t) {
						n ? n(t) : console.error(t), a.manager.itemError(e), a.manager.itemEnd(e)
					},
					o = new FileLoader(this.manager);
				o.setPath(this.path), o.setResponseType("arraybuffer"), o.setRequestHeader(this.requestHeader), o.setWithCredentials(this.withCredentials), o.load(e, (function(i) {
					try {
						a.parse(i, r, (function(i) {
							t(i), a.manager.itemEnd(e)
						}), s)
					} catch (e) {
						s(e)
					}
				}), i, s)
			},
			setDRACOLoader: function(e) {
				return this.dracoLoader = e, this
			},
			setDDSLoader: function() {
				throw new Error('THREE.GLTFLoader: "MSFT_texture_dds" no longer supported. Please update to "KHR_texture_basisu".')
			},
			setKTX2Loader: function(e) {
				return this.ktx2Loader = e, this
			},
			setMeshoptDecoder: function(e) {
				return this.meshoptDecoder = e, this
			},
			register: function(e) {
				return -1 === this.pluginCallbacks.indexOf(e) && this.pluginCallbacks.push(e), this
			},
			unregister: function(e) {
				return -1 !== this.pluginCallbacks.indexOf(e) && this.pluginCallbacks.splice(this.pluginCallbacks.indexOf(e), 1), this
			},
			parse: function(e, t, n, a) {
				var s, o = {},
					l = {};
				if ("string" == typeof e) s = e;
				else if (LoaderUtils.decodeText(new Uint8Array(e, 0, 4)) === h) {
					try {
						o[i.KHR_BINARY_GLTF] = new p(e)
					} catch (e) {
						return void(a && a(e))
					}
					s = o[i.KHR_BINARY_GLTF].content
				} else s = LoaderUtils.decodeText(new Uint8Array(e));
				var c = JSON.parse(s);
				if (void 0 === c.asset || c.asset.version[0] < 2) a && a(new Error("THREE.GLTFLoader: Unsupported asset. glTF versions >=2.0 are supported."));
				else {
					var u = new V(c, {
						path: t || this.resourcePath || "",
						crossOrigin: this.crossOrigin,
						requestHeader: this.requestHeader,
						manager: this.manager,
						ktx2Loader: this.ktx2Loader,
						meshoptDecoder: this.meshoptDecoder
					});
					u.fileLoader.setRequestHeader(this.requestHeader);
					for (var d = 0; d < this.pluginCallbacks.length; d++) {
						var g = this.pluginCallbacks[d](u);
						l[g.name] = g, o[g.name] = !0
					}
					if (c.extensionsUsed)
						for (d = 0; d < c.extensionsUsed.length; ++d) {
							var y = c.extensionsUsed[d],
								E = c.extensionsRequired || [];
							switch (y) {
								case i.KHR_MATERIALS_UNLIT:
									o[y] = new r;
									break;
								case i.KHR_MATERIALS_PBR_SPECULAR_GLOSSINESS:
									o[y] = new f;
									break;
								case i.KHR_DRACO_MESH_COMPRESSION:
									o[y] = new m(c, this.dracoLoader);
									break;
								case i.KHR_TEXTURE_TRANSFORM:
									o[y] = new A;
									break;
								case i.KHR_MESH_QUANTIZATION:
									o[y] = new v;
									break;
								default:
									E.indexOf(y) >= 0 && void 0 === l[y] && console.warn('THREE.GLTFLoader: Unknown extension "' + y + '".')
							}
						}
					u.setExtensions(o), u.setPlugins(l), u.parse(n, a)
				}
			}
		});
		var i = {
			KHR_BINARY_GLTF: "KHR_binary_glTF",
			KHR_DRACO_MESH_COMPRESSION: "KHR_draco_mesh_compression",
			KHR_LIGHTS_PUNCTUAL: "KHR_lights_punctual",
			KHR_MATERIALS_CLEARCOAT: "KHR_materials_clearcoat",
			KHR_MATERIALS_PBR_SPECULAR_GLOSSINESS: "KHR_materials_pbrSpecularGlossiness",
			KHR_MATERIALS_TRANSMISSION: "KHR_materials_transmission",
			KHR_MATERIALS_UNLIT: "KHR_materials_unlit",
			KHR_TEXTURE_BASISU: "KHR_texture_basisu",
			KHR_TEXTURE_TRANSFORM: "KHR_texture_transform",
			KHR_MESH_QUANTIZATION: "KHR_mesh_quantization",
			EXT_TEXTURE_WEBP: "EXT_texture_webp",
			EXT_MESHOPT_COMPRESSION: "EXT_meshopt_compression"
		};

		function n(e) {
			this.parser = e, this.name = i.KHR_LIGHTS_PUNCTUAL, this.cache = {
				refs: {},
				uses: {}
			}
		}

		function r() {
			this.name = i.KHR_MATERIALS_UNLIT
		}

		function a(e) {
			this.parser = e, this.name = i.KHR_MATERIALS_CLEARCOAT
		}

		function s(e) {
			this.parser = e, this.name = i.KHR_MATERIALS_TRANSMISSION
		}

		function o(e) {
			this.parser = e, this.name = i.KHR_TEXTURE_BASISU
		}

		function l(e) {
			this.parser = e, this.name = i.EXT_TEXTURE_WEBP, this.isSupported = null
		}

		function c(e) {
			this.name = i.EXT_MESHOPT_COMPRESSION, this.parser = e
		}
		n.prototype._markDefs = function() {
			for (var e = this.parser, t = this.parser.json.nodes || [], i = 0, n = t.length; i < n; i++) {
				var r = t[i];
				r.extensions && r.extensions[this.name] && void 0 !== r.extensions[this.name].light && e._addNodeRef(this.cache, r.extensions[this.name].light)
			}
		}, n.prototype._loadLight = function(e) {
			var t = this.parser,
				i = "light:" + e,
				n = t.cache.get(i);
			if (n) return n;
			var r, a = t.json,
				s = ((a.extensions && a.extensions[this.name] || {}).lights || [])[e],
				o = new Color(16777215);
			void 0 !== s.color && o.fromArray(s.color);
			var l = void 0 !== s.range ? s.range : 0;
			switch (s.type) {
				case "directional":
					(r = new DirectionalLight(o)).target.position.set(0, 0, -1), r.add(r.target);
					break;
				case "point":
					(r = new PointLight(o)).distance = l;
					break;
				case "spot":
					(r = new SpotLight(o)).distance = l, s.spot = s.spot || {}, s.spot.innerConeAngle = void 0 !== s.spot.innerConeAngle ? s.spot.innerConeAngle : 0, s.spot.outerConeAngle = void 0 !== s.spot.outerConeAngle ? s.spot.outerConeAngle : Math.PI / 4, r.angle = s.spot.outerConeAngle, r.penumbra = 1 - s.spot.innerConeAngle / s.spot.outerConeAngle, r.target.position.set(0, 0, -1), r.add(r.target);
					break;
				default:
					throw new Error("THREE.GLTFLoader: Unexpected light type: " + s.type)
			}
			return r.position.set(0, 0, 0), r.decay = 2, void 0 !== s.intensity && (r.intensity = s.intensity), r.name = t.createUniqueName(s.name || "light_" + e), n = Promise.resolve(r), t.cache.add(i, n), n
		}, n.prototype.createNodeAttachment = function(e) {
			var t = this,
				i = this.parser,
				n = i.json.nodes[e],
				r = (n.extensions && n.extensions[this.name] || {}).light;
			return void 0 === r ? null : this._loadLight(r).then((function(e) {
				return i._getNodeRef(t.cache, r, e)
			}))
		}, r.prototype.getMaterialType = function() {
			return MeshBasicMaterial
		}, r.prototype.extendParams = function(e, t, i) {
			var n = [];
			e.color = new Color(1, 1, 1), e.opacity = 1;
			var r = t.pbrMetallicRoughness;
			if (r) {
				if (Array.isArray(r.baseColorFactor)) {
					var a = r.baseColorFactor;
					e.color.fromArray(a), e.opacity = a[3]
				}
				void 0 !== r.baseColorTexture && n.push(i.assignTexture(e, "map", r.baseColorTexture))
			}
			return Promise.all(n)
		}, a.prototype.getMaterialType = function(e) {
			var t = this.parser.json.materials[e];
			return t.extensions && t.extensions[this.name] ? MeshPhysicalMaterial : null
		}, a.prototype.extendMaterialParams = function(e, t) {
			var i = this.parser,
				n = i.json.materials[e];
			if (!n.extensions || !n.extensions[this.name]) return Promise.resolve();
			var r = [],
				a = n.extensions[this.name];
			if (void 0 !== a.clearcoatFactor && (t.clearcoat = a.clearcoatFactor), void 0 !== a.clearcoatTexture && r.push(i.assignTexture(t, "clearcoatMap", a.clearcoatTexture)), void 0 !== a.clearcoatRoughnessFactor && (t.clearcoatRoughness = a.clearcoatRoughnessFactor), void 0 !== a.clearcoatRoughnessTexture && r.push(i.assignTexture(t, "clearcoatRoughnessMap", a.clearcoatRoughnessTexture)), void 0 !== a.clearcoatNormalTexture && (r.push(i.assignTexture(t, "clearcoatNormalMap", a.clearcoatNormalTexture)), void 0 !== a.clearcoatNormalTexture.scale)) {
				var s = a.clearcoatNormalTexture.scale;
				t.clearcoatNormalScale = new Vector2(s, -s)
			}
			return Promise.all(r)
		}, s.prototype.getMaterialType = function(e) {
			var t = this.parser.json.materials[e];
			return t.extensions && t.extensions[this.name] ? MeshPhysicalMaterial : null
		}, s.prototype.extendMaterialParams = function(e, t) {
			var i = this.parser,
				n = i.json.materials[e];
			if (!n.extensions || !n.extensions[this.name]) return Promise.resolve();
			var r = [],
				a = n.extensions[this.name];
			return void 0 !== a.transmissionFactor && (t.transmission = a.transmissionFactor), void 0 !== a.transmissionTexture && r.push(i.assignTexture(t, "transmissionMap", a.transmissionTexture)), Promise.all(r)
		}, o.prototype.loadTexture = function(e) {
			var t = this.parser,
				i = t.json,
				n = i.textures[e];
			if (!n.extensions || !n.extensions[this.name]) return null;
			var r = n.extensions[this.name],
				a = i.images[r.source],
				s = t.options.ktx2Loader;
			if (!s) {
				if (i.extensionsRequired && i.extensionsRequired.indexOf(this.name) >= 0) throw new Error("THREE.GLTFLoader: setKTX2Loader must be called before loading KTX2 textures");
				return null
			}
			return t.loadTextureImage(e, a, s)
		}, l.prototype.loadTexture = function(e) {
			var t = this.name,
				i = this.parser,
				n = i.json,
				r = n.textures[e];
			if (!r.extensions || !r.extensions[t]) return null;
			var a = r.extensions[t],
				s = n.images[a.source],
				o = i.textureLoader;
			if (s.uri) {
				var l = i.options.manager.getHandler(s.uri);
				null !== l && (o = l)
			}
			return this.detectSupport().then((function(r) {
				if (r) return i.loadTextureImage(e, s, o);
				if (n.extensionsRequired && n.extensionsRequired.indexOf(t) >= 0) throw new Error("THREE.GLTFLoader: WebP required by asset but unsupported.");
				return i.loadTexture(e)
			}))
		}, l.prototype.detectSupport = function() {
			return this.isSupported || (this.isSupported = new Promise((function(e) {
				var t = new Image;
				t.src = "data:image/webp;base64,UklGRiIAAABXRUJQVlA4IBYAAAAwAQCdASoBAAEADsD+JaQAA3AAAAAA", t.onload = t.onerror = function() {
					e(1 === t.height)
				}
			}))), this.isSupported
		}, c.prototype.loadBufferView = function(e) {
			var t = this.parser.json,
				i = t.bufferViews[e];
			if (i.extensions && i.extensions[this.name]) {
				var n = i.extensions[this.name],
					r = this.parser.getDependency("buffer", n.buffer),
					a = this.parser.options.meshoptDecoder;
				if (!a || !a.supported) {
					if (t.extensionsRequired && t.extensionsRequired.indexOf(this.name) >= 0) throw new Error("THREE.GLTFLoader: setMeshoptDecoder must be called before loading compressed files");
					return null
				}
				return Promise.all([r, a.ready]).then((function(e) {
					var t = n.byteOffset || 0,
						i = n.byteLength || 0,
						r = n.count,
						s = n.byteStride,
						o = new ArrayBuffer(r * s),
						l = new Uint8Array(e[0], t, i);
					return a.decodeGltfBuffer(new Uint8Array(o), r, s, l, n.mode, n.filter), o
				}))
			}
			return null
		};
		var h = "glTF",
			u = 1313821514,
			d = 5130562;

		function p(e) {
			this.name = i.KHR_BINARY_GLTF, this.content = null, this.body = null;
			var t = new DataView(e, 0, 12);
			if (this.header = {
					magic: LoaderUtils.decodeText(new Uint8Array(e.slice(0, 4))),
					version: t.getUint32(4, !0),
					length: t.getUint32(8, !0)
				}, this.header.magic !== h) throw new Error("THREE.GLTFLoader: Unsupported glTF-Binary header.");
			if (this.header.version < 2) throw new Error("THREE.GLTFLoader: Legacy binary file detected.");
			for (var n = this.header.length - 12, r = new DataView(e, 12), a = 0; a < n;) {
				var s = r.getUint32(a, !0);
				a += 4;
				var o = r.getUint32(a, !0);
				if (a += 4, o === u) {
					var l = new Uint8Array(e, 12 + a, s);
					this.content = LoaderUtils.decodeText(l)
				} else if (o === d) {
					var c = 12 + a;
					this.body = e.slice(c, c + s)
				}
				a += s
			}
			if (null === this.content) throw new Error("THREE.GLTFLoader: JSON content not found.")
		}

		function m(e, t) {
			if (!t) throw new Error("THREE.GLTFLoader: No DRACOLoader instance provided.");
			this.name = i.KHR_DRACO_MESH_COMPRESSION, this.json = e, this.dracoLoader = t, this.dracoLoader.preload()
		}

		function A() {
			this.name = i.KHR_TEXTURE_TRANSFORM
		}

		function g(e) {
			MeshStandardMaterial.call(this), this.isGLTFSpecularGlossinessMaterial = !0;
			var t = ["#ifdef USE_SPECULARMAP", "\tuniform sampler2D specularMap;", "#endif"].join("\n"),
				i = ["#ifdef USE_GLOSSINESSMAP", "\tuniform sampler2D glossinessMap;", "#endif"].join("\n"),
				n = ["vec3 specularFactor = specular;", "#ifdef USE_SPECULARMAP", "\tvec4 texelSpecular = texture2D( specularMap, vUv );", "\ttexelSpecular = sRGBToLinear( texelSpecular );", "\t// reads channel RGB, compatible with a glTF Specular-Glossiness (RGBA) texture", "\tspecularFactor *= texelSpecular.rgb;", "#endif"].join("\n"),
				r = ["float glossinessFactor = glossiness;", "#ifdef USE_GLOSSINESSMAP", "\tvec4 texelGlossiness = texture2D( glossinessMap, vUv );", "\t// reads channel A, compatible with a glTF Specular-Glossiness (RGBA) texture", "\tglossinessFactor *= texelGlossiness.a;", "#endif"].join("\n"),
				a = ["PhysicalMaterial material;", "material.diffuseColor = diffuseColor.rgb * ( 1. - max( specularFactor.r, max( specularFactor.g, specularFactor.b ) ) );", "vec3 dxy = max( abs( dFdx( geometryNormal ) ), abs( dFdy( geometryNormal ) ) );", "float geometryRoughness = max( max( dxy.x, dxy.y ), dxy.z );", "material.specularRoughness = max( 1.0 - glossinessFactor, 0.0525 ); // 0.0525 corresponds to the base mip of a 256 cubemap.", "material.specularRoughness += geometryRoughness;", "material.specularRoughness = min( material.specularRoughness, 1.0 );", "material.specularColor = specularFactor;"].join("\n"),
				s = {
					specular: {
						value: (new Color).setHex(16777215)
					},
					glossiness: {
						value: 1
					},
					specularMap: {
						value: null
					},
					glossinessMap: {
						value: null
					}
				};
			this._extraUniforms = s, this.onBeforeCompile = function(e) {
				for (var o in s) e.uniforms[o] = s[o];
				e.fragmentShader = e.fragmentShader.replace("uniform float roughness;", "uniform vec3 specular;").replace("uniform float metalness;", "uniform float glossiness;").replace("#include <roughnessmap_pars_fragment>", t).replace("#include <metalnessmap_pars_fragment>", i).replace("#include <roughnessmap_fragment>", n).replace("#include <metalnessmap_fragment>", r).replace("#include <lights_physical_fragment>", a)
			}, Object.defineProperties(this, {
				specular: {
					get: function() {
						return s.specular.value
					},
					set: function(e) {
						s.specular.value = e
					}
				},
				specularMap: {
					get: function() {
						return s.specularMap.value
					},
					set: function(e) {
						s.specularMap.value = e, e ? this.defines.USE_SPECULARMAP = "" : delete this.defines.USE_SPECULARMAP
					}
				},
				glossiness: {
					get: function() {
						return s.glossiness.value
					},
					set: function(e) {
						s.glossiness.value = e
					}
				},
				glossinessMap: {
					get: function() {
						return s.glossinessMap.value
					},
					set: function(e) {
						s.glossinessMap.value = e, e ? (this.defines.USE_GLOSSINESSMAP = "", this.defines.USE_UV = "") : (delete this.defines.USE_GLOSSINESSMAP, delete this.defines.USE_UV)
					}
				}
			}), delete this.metalness, delete this.roughness, delete this.metalnessMap, delete this.roughnessMap, this.setValues(e)
		}

		function f() {
			return {
				name: i.KHR_MATERIALS_PBR_SPECULAR_GLOSSINESS,
				specularGlossinessParams: ["color", "map", "lightMap", "lightMapIntensity", "aoMap", "aoMapIntensity", "emissive", "emissiveIntensity", "emissiveMap", "bumpMap", "bumpScale", "normalMap", "normalMapType", "displacementMap", "displacementScale", "displacementBias", "specularMap", "specular", "glossinessMap", "glossiness", "alphaMap", "envMap", "envMapIntensity", "refractionRatio"],
				getMaterialType: function() {
					return g
				},
				extendParams: function(e, t, i) {
					var n = t.extensions[this.name];
					e.color = new Color(1, 1, 1), e.opacity = 1;
					var r = [];
					if (Array.isArray(n.diffuseFactor)) {
						var a = n.diffuseFactor;
						e.color.fromArray(a), e.opacity = a[3]
					}
					if (void 0 !== n.diffuseTexture && r.push(i.assignTexture(e, "map", n.diffuseTexture)), e.emissive = new Color(0, 0, 0), e.glossiness = void 0 !== n.glossinessFactor ? n.glossinessFactor : 1, e.specular = new Color(1, 1, 1), Array.isArray(n.specularFactor) && e.specular.fromArray(n.specularFactor), void 0 !== n.specularGlossinessTexture) {
						var s = n.specularGlossinessTexture;
						r.push(i.assignTexture(e, "glossinessMap", s)), r.push(i.assignTexture(e, "specularMap", s))
					}
					return Promise.all(r)
				},
				createMaterial: function(e) {
					var t = new g(e);
					return t.fog = !0, t.color = e.color, t.map = void 0 === e.map ? null : e.map, t.lightMap = null, t.lightMapIntensity = 1, t.aoMap = void 0 === e.aoMap ? null : e.aoMap, t.aoMapIntensity = 1, t.emissive = e.emissive, t.emissiveIntensity = 1, t.emissiveMap = void 0 === e.emissiveMap ? null : e.emissiveMap, t.bumpMap = void 0 === e.bumpMap ? null : e.bumpMap, t.bumpScale = 1, t.normalMap = void 0 === e.normalMap ? null : e.normalMap, t.normalMapType = 0, e.normalScale && (t.normalScale = e.normalScale), t.displacementMap = null, t.displacementScale = 1, t.displacementBias = 0, t.specularMap = void 0 === e.specularMap ? null : e.specularMap, t.specular = e.specular, t.glossinessMap = void 0 === e.glossinessMap ? null : e.glossinessMap, t.glossiness = e.glossiness, t.alphaMap = null, t.envMap = void 0 === e.envMap ? null : e.envMap, t.envMapIntensity = 1, t.refractionRatio = .98, t
				}
			}
		}

		function v() {
			this.name = i.KHR_MESH_QUANTIZATION
		}

		function y(e, t, i, n) {
			Interpolant.call(this, e, t, i, n)
		}
		m.prototype.decodePrimitive = function(e, t) {
			var i = this.json,
				n = this.dracoLoader,
				r = e.extensions[this.name].bufferView,
				a = e.extensions[this.name].attributes,
				s = {},
				o = {},
				l = {};
			for (var c in a) {
				var h = L[c] || c.toLowerCase();
				s[h] = a[c]
			}
			for (c in e.attributes) {
				h = L[c] || c.toLowerCase();
				if (void 0 !== a[c]) {
					var u = i.accessors[e.attributes[c]],
						d = I[u.componentType];
					l[h] = d, o[h] = !0 === u.normalized
				}
			}
			return t.getDependency("bufferView", r).then((function(e) {
				return new Promise((function(t) {
					n.decodeDracoFile(e, (function(e) {
						for (var i in e.attributes) {
							var n = e.attributes[i],
								r = o[i];
							void 0 !== r && (n.normalized = r)
						}
						t(e)
					}), s, l)
				}))
			}))
		}, A.prototype.extendTexture = function(e, t) {
			return e = e.clone(), void 0 !== t.offset && e.offset.fromArray(t.offset), void 0 !== t.rotation && (e.rotation = t.rotation), void 0 !== t.scale && e.repeat.fromArray(t.scale), void 0 !== t.texCoord && console.warn('THREE.GLTFLoader: Custom UV sets in "' + this.name + '" extension not yet supported.'), e.needsUpdate = !0, e
		}, g.prototype = Object.create(MeshStandardMaterial.prototype), g.prototype.constructor = g, g.prototype.copy = function(e) {
			return MeshStandardMaterial.prototype.copy.call(this, e), this.specularMap = e.specularMap, this.specular.copy(e.specular), this.glossinessMap = e.glossinessMap, this.glossiness = e.glossiness, delete this.metalness, delete this.roughness, delete this.metalnessMap, delete this.roughnessMap, this
		}, y.prototype = Object.create(Interpolant.prototype), y.prototype.constructor = y, y.prototype.copySampleValue_ = function(e) {
			for (var t = this.resultBuffer, i = this.sampleValues, n = this.valueSize, r = e * n * 3 + n, a = 0; a !== n; a++) t[a] = i[r + a];
			return t
		}, y.prototype.beforeStart_ = y.prototype.copySampleValue_, y.prototype.afterEnd_ = y.prototype.copySampleValue_, y.prototype.interpolate_ = function(e, t, i, n) {
			for (var r = this.resultBuffer, a = this.sampleValues, s = this.valueSize, o = 2 * s, l = 3 * s, c = n - t, h = (i - t) / c, u = h * h, d = u * h, p = e * l, m = p - l, A = -2 * d + 3 * u, g = d - u, f = 1 - A, v = g - u + h, y = 0; y !== s; y++) {
				var E = a[m + y + s],
					_ = a[m + y + o] * c,
					b = a[p + y + s],
					x = a[p + y] * c;
				r[y] = f * E + v * _ + A * b + g * x
			}
			return r
		};
		var E = 0,
			_ = 1,
			b = 2,
			x = 3,
			w = 4,
			C = 5,
			S = 6,
			I = {
				5120: Int8Array,
				5121: Uint8Array,
				5122: Int16Array,
				5123: Uint16Array,
				5125: Uint32Array,
				5126: Float32Array
			},
			M = {
				9728: 1003,
				9729: 1006,
				9984: 1004,
				9985: 1007,
				9986: 1005,
				9987: 1008
			},
			T = {
				33071: 1001,
				33648: 1002,
				10497: 1e3
			},
			B = {
				SCALAR: 1,
				VEC2: 2,
				VEC3: 3,
				VEC4: 4,
				MAT2: 4,
				MAT3: 9,
				MAT4: 16
			},
			L = {
				POSITION: "position",
				NORMAL: "normal",
				TANGENT: "tangent",
				TEXCOORD_0: "uv",
				TEXCOORD_1: "uv2",
				COLOR_0: "color",
				WEIGHTS_0: "skinWeight",
				JOINTS_0: "skinIndex"
			},
			R = {
				scale: "scale",
				translation: "position",
				rotation: "quaternion",
				weights: "morphTargetInfluences"
			},
			D = {
				CUBICSPLINE: void 0,
				LINEAR: 2301,
				STEP: 2300
			},
			P = "OPAQUE",
			Q = "MASK",
			F = "BLEND";

		function O(e, t) {
			return "string" != typeof e || "" === e ? "" : (/^https?:\/\//i.test(t) && /^\//.test(e) && (t = t.replace(/(^https?:\/\/[^\/]+).*/i, "$1")), /^(https?:)?\/\//i.test(e) || /^data:.*,.*$/i.test(e) || /^blob:.*$/i.test(e) ? e : t + e)
		}

		function N(e) {
			return void 0 === e.DefaultMaterial && (e.DefaultMaterial = new MeshStandardMaterial({
				color: 16777215,
				emissive: 0,
				metalness: 1,
				roughness: 1,
				transparent: !1,
				depthTest: !0,
				side: 0
			})), e.DefaultMaterial
		}

		function k(e, t, i) {
			for (var n in i.extensions) void 0 === e[n] && (t.userData.gltfExtensions = t.userData.gltfExtensions || {}, t.userData.gltfExtensions[n] = i.extensions[n])
		}

		function U(e, t) {
			void 0 !== t.extras && ("object" == typeof t.extras ? Object.assign(e.userData, t.extras) : console.warn("THREE.GLTFLoader: Ignoring primitive type .extras, " + t.extras))
		}

		function G(e, t) {
			if (e.updateMorphTargets(), void 0 !== t.weights)
				for (var i = 0, n = t.weights.length; i < n; i++) e.morphTargetInfluences[i] = t.weights[i];
			if (t.extras && Array.isArray(t.extras.targetNames)) {
				var r = t.extras.targetNames;
				if (e.morphTargetInfluences.length === r.length) {
					e.morphTargetDictionary = {};
					for (i = 0, n = r.length; i < n; i++) e.morphTargetDictionary[r[i]] = i
				} else console.warn("THREE.GLTFLoader: Invalid extras.targetNames length. Ignoring names.")
			}
		}

		function H(e) {
			for (var t = "", i = Object.keys(e).sort(), n = 0, r = i.length; n < r; n++) t += i[n] + ":" + e[i[n]] + ";";
			return t
		}

		function V(e, i) {
			this.json = e || {}, this.extensions = {}, this.plugins = {}, this.options = i || {}, this.cache = new t, this.associations = new Map, this.primitiveCache = {}, this.meshCache = {
				refs: {},
				uses: {}
			}, this.cameraCache = {
				refs: {},
				uses: {}
			}, this.lightCache = {
				refs: {},
				uses: {}
			}, this.nodeNamesUsed = {}, "undefined" != typeof createImageBitmap && !1 === /Firefox/.test(navigator.userAgent) ? this.textureLoader = new ImageBitmapLoader(this.options.manager) : this.textureLoader = new TextureLoader(this.options.manager), this.textureLoader.setCrossOrigin(this.options.crossOrigin), this.textureLoader.setRequestHeader(this.options.requestHeader), this.fileLoader = new FileLoader(this.options.manager), this.fileLoader.setResponseType("arraybuffer"), "use-credentials" === this.options.crossOrigin && this.fileLoader.setWithCredentials(!0)
		}

		function $(e, t, i) {
			var n = t.attributes,
				r = [];

			function a(t, n) {
				return i.getDependency("accessor", t).then((function(t) {
					e.setAttribute(n, t)
				}))
			}
			for (var s in n) {
				var o = L[s] || s.toLowerCase();
				o in e.attributes || r.push(a(n[s], o))
			}
			if (void 0 !== t.indices && !e.index) {
				var l = i.getDependency("accessor", t.indices).then((function(t) {
					e.setIndex(t)
				}));
				r.push(l)
			}
			return U(e, t),
				function(e, t, i) {
					var n = t.attributes,
						r = new Box3;
					if (void 0 !== n.POSITION) {
						var a = (p = i.json.accessors[n.POSITION]).min,
							s = p.max;
						if (void 0 !== a && void 0 !== s) {
							r.set(new Vector3(a[0], a[1], a[2]), new Vector3(s[0], s[1], s[2]));
							var o = t.targets;
							if (void 0 !== o) {
								for (var l = new Vector3, c = new Vector3, h = 0, u = o.length; h < u; h++) {
									var d = o[h];
									if (void 0 !== d.POSITION) {
										var p;
										a = (p = i.json.accessors[d.POSITION]).min, s = p.max;
										void 0 !== a && void 0 !== s ? (c.setX(Math.max(Math.abs(a[0]), Math.abs(s[0]))), c.setY(Math.max(Math.abs(a[1]), Math.abs(s[1]))), c.setZ(Math.max(Math.abs(a[2]), Math.abs(s[2]))), l.max(c)) : console.warn("THREE.GLTFLoader: Missing min/max properties for accessor POSITION.")
									}
								}
								r.expandByVector(l)
							}
							e.boundingBox = r;
							var m = new Sphere;
							r.getCenter(m.center), m.radius = r.min.distanceTo(r.max) / 2, e.boundingSphere = m
						} else console.warn("THREE.GLTFLoader: Missing min/max properties for accessor POSITION.")
					}
				}(e, t, i), Promise.all(r).then((function() {
					return void 0 !== t.targets ? function(e, t, i) {
						for (var n = !1, r = !1, a = 0, s = t.length; a < s; a++) {
							if (void 0 !== (c = t[a]).POSITION && (n = !0), void 0 !== c.NORMAL && (r = !0), n && r) break
						}
						if (!n && !r) return Promise.resolve(e);
						var o = [],
							l = [];
						for (a = 0, s = t.length; a < s; a++) {
							var c = t[a];
							if (n) {
								var h = void 0 !== c.POSITION ? i.getDependency("accessor", c.POSITION) : e.attributes.position;
								o.push(h)
							}
							if (r) {
								h = void 0 !== c.NORMAL ? i.getDependency("accessor", c.NORMAL) : e.attributes.normal;
								l.push(h)
							}
						}
						return Promise.all([Promise.all(o), Promise.all(l)]).then((function(t) {
							var i = t[0],
								a = t[1];
							return n && (e.morphAttributes.position = i), r && (e.morphAttributes.normal = a), e.morphTargetsRelative = !0, e
						}))
					}(e, t.targets, i) : e
				}))
		}

		function z(e, t) {
			var i = e.getIndex();
			if (null === i) {
				var n = [],
					r = e.getAttribute("position");
				if (void 0 === r) return console.error("THREE.GLTFLoader.toTrianglesDrawMode(): Undefined position attribute. Processing not possible."), e;
				for (var a = 0; a < r.count; a++) n.push(a);
				e.setIndex(n), i = e.getIndex()
			}
			var s = i.count - 2,
				o = [];
			if (2 === t)
				for (a = 1; a <= s; a++) o.push(i.getX(0)), o.push(i.getX(a)), o.push(i.getX(a + 1));
			else
				for (a = 0; a < s; a++) a % 2 == 0 ? (o.push(i.getX(a)), o.push(i.getX(a + 1)), o.push(i.getX(a + 2))) : (o.push(i.getX(a + 2)), o.push(i.getX(a + 1)), o.push(i.getX(a)));
			o.length / 3 !== s && console.error("THREE.GLTFLoader.toTrianglesDrawMode(): Unable to generate correct amount of triangles.");
			var l = e.clone();
			return l.setIndex(o), l
		}
		return V.prototype.setExtensions = function(e) {
			this.extensions = e
		}, V.prototype.setPlugins = function(e) {
			this.plugins = e
		}, V.prototype.parse = function(e, t) {
			var i = this,
				n = this.json,
				r = this.extensions;
			this.cache.removeAll(), this._invokeAll((function(e) {
				return e._markDefs && e._markDefs()
			})), Promise.all(this._invokeAll((function(e) {
				return e.beforeRoot && e.beforeRoot()
			}))).then((function() {
				return Promise.all([i.getDependencies("scene"), i.getDependencies("animation"), i.getDependencies("camera")])
			})).then((function(t) {
				var a = {
					scene: t[0][n.scene || 0],
					scenes: t[0],
					animations: t[1],
					cameras: t[2],
					asset: n.asset,
					parser: i,
					userData: {}
				};
				k(r, a, n), U(a, n), Promise.all(i._invokeAll((function(e) {
					return e.afterRoot && e.afterRoot(a)
				}))).then((function() {
					e(a)
				}))
			})).catch(t)
		}, V.prototype._markDefs = function() {
			for (var e = this.json.nodes || [], t = this.json.skins || [], i = this.json.meshes || [], n = 0, r = t.length; n < r; n++)
				for (var a = t[n].joints, s = 0, o = a.length; s < o; s++) e[a[s]].isBone = !0;
			for (var l = 0, c = e.length; l < c; l++) {
				var h = e[l];
				void 0 !== h.mesh && (this._addNodeRef(this.meshCache, h.mesh), void 0 !== h.skin && (i[h.mesh].isSkinnedMesh = !0)), void 0 !== h.camera && this._addNodeRef(this.cameraCache, h.camera)
			}
		}, V.prototype._addNodeRef = function(e, t) {
			void 0 !== t && (void 0 === e.refs[t] && (e.refs[t] = e.uses[t] = 0), e.refs[t]++)
		}, V.prototype._getNodeRef = function(e, t, i) {
			if (e.refs[t] <= 1) return i;
			var n = i.clone();
			return n.name += "_instance_" + e.uses[t]++, n
		}, V.prototype._invokeOne = function(e) {
			var t = Object.values(this.plugins);
			t.push(this);
			for (var i = 0; i < t.length; i++) {
				var n = e(t[i]);
				if (n) return n
			}
		}, V.prototype._invokeAll = function(e) {
			var t = Object.values(this.plugins);
			t.unshift(this);
			for (var i = [], n = 0; n < t.length; n++) {
				var r = e(t[n]);
				r && i.push(r)
			}
			return i
		}, V.prototype.getDependency = function(e, t) {
			var i = e + ":" + t,
				n = this.cache.get(i);
			if (!n) {
				switch (e) {
					case "scene":
						n = this.loadScene(t);
						break;
					case "node":
						n = this.loadNode(t);
						break;
					case "mesh":
						n = this._invokeOne((function(e) {
							return e.loadMesh && e.loadMesh(t)
						}));
						break;
					case "accessor":
						n = this.loadAccessor(t);
						break;
					case "bufferView":
						n = this._invokeOne((function(e) {
							return e.loadBufferView && e.loadBufferView(t)
						}));
						break;
					case "buffer":
						n = this.loadBuffer(t);
						break;
					case "material":
						n = this._invokeOne((function(e) {
							return e.loadMaterial && e.loadMaterial(t)
						}));
						break;
					case "texture":
						n = this._invokeOne((function(e) {
							return e.loadTexture && e.loadTexture(t)
						}));
						break;
					case "skin":
						n = this.loadSkin(t);
						break;
					case "animation":
						n = this.loadAnimation(t);
						break;
					case "camera":
						n = this.loadCamera(t);
						break;
					default:
						throw new Error("Unknown type: " + e)
				}
				this.cache.add(i, n)
			}
			return n
		}, V.prototype.getDependencies = function(e) {
			var t = this.cache.get(e);
			if (!t) {
				var i = this,
					n = this.json[e + ("mesh" === e ? "es" : "s")] || [];
				t = Promise.all(n.map((function(t, n) {
					return i.getDependency(e, n)
				}))), this.cache.add(e, t)
			}
			return t
		}, V.prototype.loadBuffer = function(e) {
			var t = this.json.buffers[e],
				n = this.fileLoader;
			if (t.type && "arraybuffer" !== t.type) throw new Error("THREE.GLTFLoader: " + t.type + " buffer type is not supported.");
			if (void 0 === t.uri && 0 === e) return Promise.resolve(this.extensions[i.KHR_BINARY_GLTF].body);
			var r = this.options;
			return new Promise((function(e, i) {
				n.load(O(t.uri, r.path), e, void 0, (function() {
					i(new Error('THREE.GLTFLoader: Failed to load buffer "' + t.uri + '".'))
				}))
			}))
		}, V.prototype.loadBufferView = function(e) {
			var t = this.json.bufferViews[e];
			return this.getDependency("buffer", t.buffer).then((function(e) {
				var i = t.byteLength || 0,
					n = t.byteOffset || 0;
				return e.slice(n, n + i)
			}))
		}, V.prototype.loadAccessor = function(e) {
			var t = this,
				i = this.json,
				n = this.json.accessors[e];
			if (void 0 === n.bufferView && void 0 === n.sparse) return Promise.resolve(null);
			var r = [];
			return void 0 !== n.bufferView ? r.push(this.getDependency("bufferView", n.bufferView)) : r.push(null), void 0 !== n.sparse && (r.push(this.getDependency("bufferView", n.sparse.indices.bufferView)), r.push(this.getDependency("bufferView", n.sparse.values.bufferView))), Promise.all(r).then((function(e) {
				var r, a = e[0],
					s = B[n.type],
					o = I[n.componentType],
					l = o.BYTES_PER_ELEMENT,
					c = l * s,
					h = n.byteOffset || 0,
					u = void 0 !== n.bufferView ? i.bufferViews[n.bufferView].byteStride : void 0,
					d = !0 === n.normalized;
				if (u && u !== c) {
					var p = Math.floor(h / u),
						m = "InterleavedBuffer:" + n.bufferView + ":" + n.componentType + ":" + p + ":" + n.count,
						A = t.cache.get(m);
					A || (A = new InterleavedBuffer(new o(a, p * u, n.count * u / l), u / l), t.cache.add(m, A)), r = new InterleavedBufferAttribute(A, s, h % u / l, d)
				} else r = new BufferAttribute(null === a ? new o(n.count * s) : new o(a, h, n.count * s), s, d);
				if (void 0 !== n.sparse) {
					var g = B.SCALAR,
						f = I[n.sparse.indices.componentType],
						v = n.sparse.indices.byteOffset || 0,
						y = n.sparse.values.byteOffset || 0,
						E = new f(e[1], v, n.sparse.count * g),
						_ = new o(e[2], y, n.sparse.count * s);
					null !== a && (r = new BufferAttribute(r.array.slice(), r.itemSize, r.normalized));
					for (var b = 0, x = E.length; b < x; b++) {
						var w = E[b];
						if (r.setX(w, _[b * s]), s >= 2 && r.setY(w, _[b * s + 1]), s >= 3 && r.setZ(w, _[b * s + 2]), s >= 4 && r.setW(w, _[b * s + 3]), s >= 5) throw new Error("THREE.GLTFLoader: Unsupported itemSize in sparse BufferAttribute.")
					}
				}
				return r
			}))
		}, V.prototype.loadTexture = function(e) {
			var t = this.json,
				i = this.options,
				n = t.textures[e],
				r = t.images[n.source],
				a = this.textureLoader;
			if (r.uri) {
				var s = i.manager.getHandler(r.uri);
				null !== s && (a = s)
			}
			return this.loadTextureImage(e, r, a)
		}, V.prototype.loadTextureImage = function(e, t, i) {
			var n = this,
				r = this.json,
				a = this.options,
				s = r.textures[e],
				o = self.URL || self.webkitURL,
				l = t.uri,
				c = !1,
				h = !0;
			if ("image/jpeg" === t.mimeType && (h = !1), void 0 !== t.bufferView) l = n.getDependency("bufferView", t.bufferView).then((function(e) {
				if ("image/png" === t.mimeType) {
					var i = new DataView(e, 25, 1).getUint8(0, !1);
					h = 6 === i || 4 === i || 3 === i
				}
				c = !0;
				var n = new Blob([e], {
					type: t.mimeType
				});
				return l = o.createObjectURL(n)
			}));
			else if (void 0 === t.uri) throw new Error("THREE.GLTFLoader: Image " + e + " is missing URI and bufferView");
			return Promise.resolve(l).then((function(e) {
				return new Promise((function(t, n) {
					var r = t;
					!0 === i.isImageBitmapLoader && (r = function(e) {
						t(new CanvasTexture(e))
					}), i.load(O(e, a.path), r, void 0, n)
				}))
			})).then((function(t) {
				!0 === c && o.revokeObjectURL(l), t.flipY = !1, s.name && (t.name = s.name), h || (t.format = 1022);
				var i = (r.samplers || {})[s.sampler] || {};
				return t.magFilter = M[i.magFilter] || 1006, t.minFilter = M[i.minFilter] || 1008, t.wrapS = T[i.wrapS] || 1e3, t.wrapT = T[i.wrapT] || 1e3, n.associations.set(t, {
					type: "textures",
					index: e
				}), t
			}))
		}, V.prototype.assignTexture = function(e, t, n) {
			var r = this;
			return this.getDependency("texture", n.index).then((function(a) {
				if (void 0 === n.texCoord || 0 == n.texCoord || "aoMap" === t && 1 == n.texCoord || console.warn("THREE.GLTFLoader: Custom UV set " + n.texCoord + " for texture " + t + " not yet supported."), r.extensions[i.KHR_TEXTURE_TRANSFORM]) {
					var s = void 0 !== n.extensions ? n.extensions[i.KHR_TEXTURE_TRANSFORM] : void 0;
					if (s) {
						var o = r.associations.get(a);
						a = r.extensions[i.KHR_TEXTURE_TRANSFORM].extendTexture(a, s), r.associations.set(a, o)
					}
				}
				e[t] = a
			}))
		}, V.prototype.assignFinalMaterial = function(e) {
			var t = e.geometry,
				i = e.material,
				n = void 0 !== t.attributes.tangent,
				r = void 0 !== t.attributes.color,
				a = void 0 === t.attributes.normal,
				s = !0 === e.isSkinnedMesh,
				o = Object.keys(t.morphAttributes).length > 0,
				l = o && void 0 !== t.morphAttributes.normal;
			if (e.isPoints) {
				var c = "PointsMaterial:" + i.uuid,
					h = this.cache.get(c);
				h || (h = new PointsMaterial, Material$1.prototype.copy.call(h, i), h.color.copy(i.color), h.map = i.map, h.sizeAttenuation = !1, this.cache.add(c, h)), i = h
			} else if (e.isLine) {
				c = "LineBasicMaterial:" + i.uuid;
				var u = this.cache.get(c);
				u || (u = new LineBasicMaterial, Material$1.prototype.copy.call(u, i), u.color.copy(i.color), this.cache.add(c, u)), i = u
			}
			if (n || r || a || s || o) {
				c = "ClonedMaterial:" + i.uuid + ":";
				i.isGLTFSpecularGlossinessMaterial && (c += "specular-glossiness:"), s && (c += "skinning:"), n && (c += "vertex-tangents:"), r && (c += "vertex-colors:"), a && (c += "flat-shading:"), o && (c += "morph-targets:"), l && (c += "morph-normals:");
				var d = this.cache.get(c);
				d || (d = i.clone(), s && (d.skinning = !0), r && (d.vertexColors = !0), a && (d.flatShading = !0), o && (d.morphTargets = !0), l && (d.morphNormals = !0), n && (d.vertexTangents = !0, d.normalScale && (d.normalScale.y *= -1), d.clearcoatNormalScale && (d.clearcoatNormalScale.y *= -1)), this.cache.add(c, d), this.associations.set(d, this.associations.get(i))), i = d
			}
			i.aoMap && void 0 === t.attributes.uv2 && void 0 !== t.attributes.uv && t.setAttribute("uv2", t.attributes.uv), e.material = i
		}, V.prototype.getMaterialType = function() {
			return MeshStandardMaterial
		}, V.prototype.loadMaterial = function(e) {
			var t, n = this,
				r = this.json,
				a = this.extensions,
				s = r.materials[e],
				o = {},
				l = s.extensions || {},
				c = [];
			if (l[i.KHR_MATERIALS_PBR_SPECULAR_GLOSSINESS]) {
				var h = a[i.KHR_MATERIALS_PBR_SPECULAR_GLOSSINESS];
				t = h.getMaterialType(), c.push(h.extendParams(o, s, n))
			} else if (l[i.KHR_MATERIALS_UNLIT]) {
				var u = a[i.KHR_MATERIALS_UNLIT];
				t = u.getMaterialType(), c.push(u.extendParams(o, s, n))
			} else {
				var d = s.pbrMetallicRoughness || {};
				if (o.color = new Color(1, 1, 1), o.opacity = 1, Array.isArray(d.baseColorFactor)) {
					var p = d.baseColorFactor;
					o.color.fromArray(p), o.opacity = p[3]
				}
				void 0 !== d.baseColorTexture && c.push(n.assignTexture(o, "map", d.baseColorTexture)), o.metalness = void 0 !== d.metallicFactor ? d.metallicFactor : 1, o.roughness = void 0 !== d.roughnessFactor ? d.roughnessFactor : 1, void 0 !== d.metallicRoughnessTexture && (c.push(n.assignTexture(o, "metalnessMap", d.metallicRoughnessTexture)), c.push(n.assignTexture(o, "roughnessMap", d.metallicRoughnessTexture))), t = this._invokeOne((function(t) {
					return t.getMaterialType && t.getMaterialType(e)
				})), c.push(Promise.all(this._invokeAll((function(t) {
					return t.extendMaterialParams && t.extendMaterialParams(e, o)
				}))))
			}!0 === s.doubleSided && (o.side = 2);
			var m = s.alphaMode || P;
			return m === F ? (o.transparent = !0, o.depthWrite = !1) : (o.transparent = !1, m === Q && (o.alphaTest = void 0 !== s.alphaCutoff ? s.alphaCutoff : .5)), void 0 !== s.normalTexture && t !== MeshBasicMaterial && (c.push(n.assignTexture(o, "normalMap", s.normalTexture)), o.normalScale = new Vector2(1, -1), void 0 !== s.normalTexture.scale && o.normalScale.set(s.normalTexture.scale, -s.normalTexture.scale)), void 0 !== s.occlusionTexture && t !== MeshBasicMaterial && (c.push(n.assignTexture(o, "aoMap", s.occlusionTexture)), void 0 !== s.occlusionTexture.strength && (o.aoMapIntensity = s.occlusionTexture.strength)), void 0 !== s.emissiveFactor && t !== MeshBasicMaterial && (o.emissive = (new Color).fromArray(s.emissiveFactor)), void 0 !== s.emissiveTexture && t !== MeshBasicMaterial && c.push(n.assignTexture(o, "emissiveMap", s.emissiveTexture)), Promise.all(c).then((function() {
				var r;
				return r = t === g ? a[i.KHR_MATERIALS_PBR_SPECULAR_GLOSSINESS].createMaterial(o) : new t(o), s.name && (r.name = s.name), r.map && (r.map.encoding = 3001), r.emissiveMap && (r.emissiveMap.encoding = 3001), U(r, s), n.associations.set(r, {
					type: "materials",
					index: e
				}), s.extensions && k(a, r, s), r
			}))
		}, V.prototype.createUniqueName = function(e) {
			for (var t = PropertyBinding.sanitizeNodeName(e || ""), i = t, n = 1; this.nodeNamesUsed[i]; ++n) i = t + "_" + n;
			return this.nodeNamesUsed[i] = !0, i
		}, V.prototype.loadGeometries = function(e) {
			var t = this,
				n = this.extensions,
				r = this.primitiveCache;

			function a(e) {
				return n[i.KHR_DRACO_MESH_COMPRESSION].decodePrimitive(e, t).then((function(i) {
					return $(i, e, t)
				}))
			}
			for (var s, o, l = [], c = 0, h = e.length; c < h; c++) {
				var u, d = e[c],
					p = (o = void 0, (o = (s = d).extensions && s.extensions[i.KHR_DRACO_MESH_COMPRESSION]) ? "draco:" + o.bufferView + ":" + o.indices + ":" + H(o.attributes) : s.indices + ":" + H(s.attributes) + ":" + s.mode),
					m = r[p];
				if (m) l.push(m.promise);
				else u = d.extensions && d.extensions[i.KHR_DRACO_MESH_COMPRESSION] ? a(d) : $(new BufferGeometry, d, t), r[p] = {
					primitive: d,
					promise: u
				}, l.push(u)
			}
			return Promise.all(l)
		}, V.prototype.loadMesh = function(e) {
			for (var t = this, i = this.json, n = this.extensions, r = i.meshes[e], a = r.primitives, s = [], o = 0, l = a.length; o < l; o++) {
				var c = void 0 === a[o].material ? N(this.cache) : this.getDependency("material", a[o].material);
				s.push(c)
			}
			return s.push(t.loadGeometries(a)), Promise.all(s).then((function(i) {
				for (var s = i.slice(0, i.length - 1), o = i[i.length - 1], l = [], c = 0, h = o.length; c < h; c++) {
					var u, d = o[c],
						p = a[c],
						m = s[c];
					if (p.mode === w || p.mode === C || p.mode === S || void 0 === p.mode) !0 !== (u = !0 === r.isSkinnedMesh ? new SkinnedMesh(d, m) : new Mesh(d, m)).isSkinnedMesh || u.geometry.attributes.skinWeight.normalized || u.normalizeSkinWeights(), p.mode === C ? u.geometry = z(u.geometry, 1) : p.mode === S && (u.geometry = z(u.geometry, 2));
					else if (p.mode === _) u = new LineSegments(d, m);
					else if (p.mode === x) u = new Line(d, m);
					else if (p.mode === b) u = new LineLoop(d, m);
					else {
						if (p.mode !== E) throw new Error("THREE.GLTFLoader: Primitive mode unsupported: " + p.mode);
						u = new Points(d, m)
					}
					Object.keys(u.geometry.morphAttributes).length > 0 && G(u, r), u.name = t.createUniqueName(r.name || "mesh_" + e), U(u, r), p.extensions && k(n, u, p), t.assignFinalMaterial(u), l.push(u)
				}
				if (1 === l.length) return l[0];
				var A = new Group;
				for (c = 0, h = l.length; c < h; c++) A.add(l[c]);
				return A
			}))
		}, V.prototype.loadCamera = function(e) {
			var t, i = this.json.cameras[e],
				n = i[i.type];
			if (n) return "perspective" === i.type ? t = new PerspectiveCamera(MathUtils.radToDeg(n.yfov), n.aspectRatio || 1, n.znear || 1, n.zfar || 2e6) : "orthographic" === i.type && (t = new OrthographicCamera(-n.xmag, n.xmag, n.ymag, -n.ymag, n.znear, n.zfar)), i.name && (t.name = this.createUniqueName(i.name)), U(t, i), Promise.resolve(t);
			console.warn("THREE.GLTFLoader: Missing camera parameters.")
		}, V.prototype.loadSkin = function(e) {
			var t = this.json.skins[e],
				i = {
					joints: t.joints
				};
			return void 0 === t.inverseBindMatrices ? Promise.resolve(i) : this.getDependency("accessor", t.inverseBindMatrices).then((function(e) {
				return i.inverseBindMatrices = e, i
			}))
		}, V.prototype.loadAnimation = function(e) {
			for (var t = this.json.animations[e], i = [], n = [], r = [], a = [], s = [], o = 0, l = t.channels.length; o < l; o++) {
				var c = t.channels[o],
					h = t.samplers[c.sampler],
					u = c.target,
					d = void 0 !== u.node ? u.node : u.id,
					p = void 0 !== t.parameters ? t.parameters[h.input] : h.input,
					m = void 0 !== t.parameters ? t.parameters[h.output] : h.output;
				i.push(this.getDependency("node", d)), n.push(this.getDependency("accessor", p)), r.push(this.getDependency("accessor", m)), a.push(h), s.push(u)
			}
			return Promise.all([Promise.all(i), Promise.all(n), Promise.all(r), Promise.all(a), Promise.all(s)]).then((function(i) {
				for (var n = i[0], r = i[1], a = i[2], s = i[3], o = i[4], l = [], c = 0, h = n.length; c < h; c++) {
					var u = n[c],
						d = r[c],
						p = a[c],
						m = s[c],
						A = o[c];
					if (void 0 !== u) {
						var g;
						switch (u.updateMatrix(), u.matrixAutoUpdate = !0, R[A.path]) {
							case R.weights:
								g = NumberKeyframeTrack;
								break;
							case R.rotation:
								g = QuaternionKeyframeTrack;
								break;
							case R.position:
							case R.scale:
							default:
								g = VectorKeyframeTrack
						}
						var f = u.name ? u.name : u.uuid,
							v = void 0 !== m.interpolation ? D[m.interpolation] : 2301,
							E = [];
						R[A.path] === R.weights ? u.traverse((function(e) {
							!0 === e.isMesh && e.morphTargetInfluences && E.push(e.name ? e.name : e.uuid)
						})) : E.push(f);
						var _ = p.array;
						if (p.normalized) {
							var b;
							if (_.constructor === Int8Array) b = 1 / 127;
							else if (_.constructor === Uint8Array) b = 1 / 255;
							else if (_.constructor == Int16Array) b = 1 / 32767;
							else {
								if (_.constructor !== Uint16Array) throw new Error("THREE.GLTFLoader: Unsupported output accessor component type.");
								b = 1 / 65535
							}
							for (var x = new Float32Array(_.length), w = 0, C = _.length; w < C; w++) x[w] = _[w] * b;
							_ = x
						}
						for (w = 0, C = E.length; w < C; w++) {
							var S = new g(E[w] + "." + R[A.path], d.array, _, v);
							"CUBICSPLINE" === m.interpolation && (S.createInterpolant = function(e) {
								return new y(this.times, this.values, this.getValueSize() / 3, e)
							}, S.createInterpolant.isInterpolantFactoryMethodGLTFCubicSpline = !0), l.push(S)
						}
					}
				}
				var I = t.name ? t.name : "animation_" + e;
				return new AnimationClip(I, void 0, l)
			}))
		}, V.prototype.loadNode = function(e) {
			var t, i = this.json,
				n = this.extensions,
				r = this,
				a = i.nodes[e],
				s = a.name ? r.createUniqueName(a.name) : "";
			return (t = [], void 0 !== a.mesh && t.push(r.getDependency("mesh", a.mesh).then((function(e) {
				var t = r._getNodeRef(r.meshCache, a.mesh, e);
				return void 0 !== a.weights && t.traverse((function(e) {
					if (e.isMesh)
						for (var t = 0, i = a.weights.length; t < i; t++) e.morphTargetInfluences[t] = a.weights[t]
				})), t
			}))), void 0 !== a.camera && t.push(r.getDependency("camera", a.camera).then((function(e) {
				return r._getNodeRef(r.cameraCache, a.camera, e)
			}))), r._invokeAll((function(t) {
				return t.createNodeAttachment && t.createNodeAttachment(e)
			})).forEach((function(e) {
				t.push(e)
			})), Promise.all(t)).then((function(t) {
				var i;
				if ((i = !0 === a.isBone ? new Bone : t.length > 1 ? new Group : 1 === t.length ? t[0] : new Object3D) !== t[0])
					for (var o = 0, l = t.length; o < l; o++) i.add(t[o]);
				if (a.name && (i.userData.name = a.name, i.name = s), U(i, a), a.extensions && k(n, i, a), void 0 !== a.matrix) {
					var c = new Matrix4;
					c.fromArray(a.matrix), i.applyMatrix4(c)
				} else void 0 !== a.translation && i.position.fromArray(a.translation), void 0 !== a.rotation && i.quaternion.fromArray(a.rotation), void 0 !== a.scale && i.scale.fromArray(a.scale);
				return r.associations.set(i, {
					type: "nodes",
					index: e
				}), i
			}))
		}, V.prototype.loadScene = function() {
			function e(t, i, n, r) {
				var a = n.nodes[t];
				return r.getDependency("node", t).then((function(e) {
					return void 0 === a.skin ? e : r.getDependency("skin", a.skin).then((function(e) {
						for (var i = [], n = 0, a = (t = e).joints.length; n < a; n++) i.push(r.getDependency("node", t.joints[n]));
						return Promise.all(i)
					})).then((function(i) {
						return e.traverse((function(e) {
							if (e.isMesh) {
								for (var n = [], r = [], a = 0, s = i.length; a < s; a++) {
									var o = i[a];
									if (o) {
										n.push(o);
										var l = new Matrix4;
										void 0 !== t.inverseBindMatrices && l.fromArray(t.inverseBindMatrices.array, 16 * a), r.push(l)
									} else console.warn('THREE.GLTFLoader: Joint "%s" could not be found.', t.joints[a])
								}
								e.bind(new Skeleton(n, r), e.matrixWorld)
							}
						})), e
					}));
					var t
				})).then((function(t) {
					i.add(t);
					var s = [];
					if (a.children)
						for (var o = a.children, l = 0, c = o.length; l < c; l++) {
							var h = o[l];
							s.push(e(h, t, n, r))
						}
					return Promise.all(s)
				}))
			}
			return function(t) {
				var i = this.json,
					n = this.extensions,
					r = this.json.scenes[t],
					a = new Group;
				r.name && (a.name = this.createUniqueName(r.name)), U(a, r), r.extensions && k(n, a, r);
				for (var s = r.nodes || [], o = [], l = 0, c = s.length; l < c; l++) o.push(e(s[l], a, i, this));
				return Promise.all(o).then((function() {
					return a
				}))
			}
		}(), e
	}(),
	BasisTextureLoader = function(e) {
		Loader.call(this, e), this.transcoderPath = "", this.transcoderBinary = null, this.transcoderPending = null, this.workerLimit = 4, this.workerPool = [], this.workerNextTaskID = 1, this.workerSourceURL = "", this.workerConfig = null
	};
let init, instance, heap;
BasisTextureLoader.taskCache = new WeakMap, BasisTextureLoader.prototype = Object.assign(Object.create(Loader.prototype), {
	constructor: BasisTextureLoader,
	setTranscoderPath: function(e) {
		return this.transcoderPath = e, this
	},
	setWorkerLimit: function(e) {
		return this.workerLimit = e, this
	},
	detectSupport: function(e) {
		return this.workerConfig = {
			astcSupported: e.extensions.has("WEBGL_compressed_texture_astc"),
			etc1Supported: e.extensions.has("WEBGL_compressed_texture_etc1"),
			etc2Supported: e.extensions.has("WEBGL_compressed_texture_etc"),
			dxtSupported: e.extensions.has("WEBGL_compressed_texture_s3tc"),
			bptcSupported: e.extensions.has("EXT_texture_compression_bptc"),
			pvrtcSupported: e.extensions.has("WEBGL_compressed_texture_pvrtc") || e.extensions.has("WEBKIT_WEBGL_compressed_texture_pvrtc")
		}, this
	},
	load: function(e, t, i, n) {
		var r = new FileLoader(this.manager);
		r.setResponseType("arraybuffer"), r.setWithCredentials(this.withCredentials);
		var a = new CompressedTexture;
		return r.load(e, e => {
			if (BasisTextureLoader.taskCache.has(e)) return BasisTextureLoader.taskCache.get(e).promise.then(t).catch(n);
			this._createTexture([e]).then((function(e) {
				a.copy(e), a.needsUpdate = !0, t && t(a)
			})).catch(n)
		}, i, n), a
	},
	parseInternalAsync: function(e) {
		for (var {
				levels: t
			} = e, i = new Set, n = 0; n < t.length; n++) i.add(t[n].data.buffer);
		return this._createTexture(Array.from(i), {
			...e,
			lowLevel: !0
		})
	},
	_createTexture: function(e, t) {
		for (var i, n, r = t || {}, a = 0, s = 0; s < e.length; s++) a += e[s].byteLength;
		var o = this._allocateWorker(a).then(t => (i = t, n = this.workerNextTaskID++, new Promise((t, a) => {
			i._callbacks[n] = {
				resolve: t,
				reject: a
			}, i.postMessage({
				type: "transcode",
				id: n,
				buffers: e,
				taskConfig: r
			}, e)
		}))).then(e => {
			var {
				mipmaps: t,
				width: i,
				height: n,
				format: r
			} = e, a = new CompressedTexture(t, i, n, r, 1009);
			return a.minFilter = 1 === t.length ? 1006 : 1008, a.magFilter = 1006, a.generateMipmaps = !1, a.needsUpdate = !0, a
		});
		return o.catch(() => !0).then(() => {
			i && n && (i._taskLoad -= a, delete i._callbacks[n])
		}), BasisTextureLoader.taskCache.set(e[0], {
			promise: o
		}), o
	},
	_initTranscoder: function() {
		if (!this.transcoderPending) {
			var e = new FileLoader(this.manager);
			e.setPath(this.transcoderPath), e.setWithCredentials(this.withCredentials);
			var t = new Promise((t, i) => {
					e.load("basis_transcoder.js", t, void 0, i)
				}),
				i = new FileLoader(this.manager);
			i.setPath(this.transcoderPath), i.setResponseType("arraybuffer"), i.setWithCredentials(this.withCredentials);
			var n = new Promise((e, t) => {
				i.load("basis_transcoder.wasm", e, void 0, t)
			});
			this.transcoderPending = Promise.all([t, n]).then(([e, t]) => {
				var i = BasisTextureLoader.BasisWorker.toString(),
					n = ["/* constants */", "var _EngineFormat = " + JSON.stringify(BasisTextureLoader.EngineFormat), "var _TranscoderFormat = " + JSON.stringify(BasisTextureLoader.TranscoderFormat), "var _BasisFormat = " + JSON.stringify(BasisTextureLoader.BasisFormat), "/* basis_transcoder.js */", e, "/* worker */", i.substring(i.indexOf("{") + 1, i.lastIndexOf("}"))].join("\n");
				this.workerSourceURL = URL.createObjectURL(new Blob([n])), this.transcoderBinary = t
			})
		}
		return this.transcoderPending
	},
	_allocateWorker: function(e) {
		return this._initTranscoder().then(() => {
			var t;
			this.workerPool.length < this.workerLimit ? ((t = new Worker(this.workerSourceURL))._callbacks = {}, t._taskLoad = 0, t.postMessage({
				type: "init",
				config: this.workerConfig,
				transcoderBinary: this.transcoderBinary
			}), t.onmessage = function(e) {
				var i = e.data;
				switch (i.type) {
					case "transcode":
						t._callbacks[i.id].resolve(i);
						break;
					case "error":
						t._callbacks[i.id].reject(i);
						break;
					default:
						console.error('THREE.BasisTextureLoader: Unexpected message, "' + i.type + '"')
				}
			}, this.workerPool.push(t)) : this.workerPool.sort((function(e, t) {
				return e._taskLoad > t._taskLoad ? -1 : 1
			}));
			return (t = this.workerPool[this.workerPool.length - 1])._taskLoad += e, t
		})
	},
	dispose: function() {
		for (var e = 0; e < this.workerPool.length; e++) this.workerPool[e].terminate();
		return this.workerPool.length = 0, this
	}
}), BasisTextureLoader.BasisFormat = {
	ETC1S: 0,
	UASTC_4x4: 1
}, BasisTextureLoader.TranscoderFormat = {
	ETC1: 0,
	ETC2: 1,
	BC1: 2,
	BC3: 3,
	BC4: 4,
	BC5: 5,
	BC7_M6_OPAQUE_ONLY: 6,
	BC7_M5: 7,
	PVRTC1_4_RGB: 8,
	PVRTC1_4_RGBA: 9,
	ASTC_4x4: 10,
	ATC_RGB: 11,
	ATC_RGBA_INTERPOLATED_ALPHA: 12,
	RGBA32: 13,
	RGB565: 14,
	BGR565: 15,
	RGBA4444: 16
}, BasisTextureLoader.EngineFormat = {
	RGBAFormat: 1023,
	RGBA_ASTC_4x4_Format: 37808,
	RGBA_BPTC_Format: 36492,
	RGBA_ETC2_EAC_Format: 37496,
	RGBA_PVRTC_4BPPV1_Format: 35842,
	RGBA_S3TC_DXT5_Format: 33779,
	RGB_ETC1_Format: 36196,
	RGB_ETC2_Format: 37492,
	RGB_PVRTC_4BPPV1_Format: 35840,
	RGB_S3TC_DXT1_Format: 33776
}, BasisTextureLoader.BasisWorker = function() {
	var e, t, i, n = _EngineFormat,
		r = _TranscoderFormat,
		a = _BasisFormat;
	onmessage = function(n) {
		var r, s = n.data;
		switch (s.type) {
			case "init":
				e = s.config, r = s.transcoderBinary, t = new Promise(e => {
					i = {
						wasmBinary: r,
						onRuntimeInitialized: e
					}, BASIS(i)
				}).then(() => {
					i.initializeBasis()
				});
				break;
			case "transcode":
				t.then(() => {
					try {
						for (var {
								width: e,
								height: t,
								hasAlpha: n,
								mipmaps: r,
								format: o
							} = s.taskConfig.lowLevel ? function(e) {
								var {
									basisFormat: t,
									width: n,
									height: r,
									hasAlpha: s
								} = e, {
									transcoderFormat: o,
									engineFormat: l
								} = c(t, n, r, s), m = i.getBytesPerBlockOrPixel(o);
								h(i.isFormatSupported(o), "THREE.BasisTextureLoader: Unsupported format.");
								var A = [];
								if (t === a.ETC1S) {
									var g = new i.LowLevelETC1SImageTranscoder,
										{
											endpointCount: f,
											endpointsData: v,
											selectorCount: y,
											selectorsData: E,
											tablesData: _
										} = e.globalData;
									try {
										h(g.decodePalettes(f, v, y, E), "THREE.BasisTextureLoader: decodePalettes() failed."), h(g.decodeTables(_), "THREE.BasisTextureLoader: decodeTables() failed.");
										for (var b = 0; b < e.levels.length; b++) {
											var x = e.levels[b],
												w = e.globalData.imageDescs[b],
												C = p(o, x.width, x.height),
												S = new Uint8Array(C);
											h(g.transcodeImage(o, S, C / m, x.data, u(o, x.width), d(o, x.height), x.width, x.height, x.index, w.rgbSliceByteOffset, w.rgbSliceByteLength, w.alphaSliceByteOffset, w.alphaSliceByteLength, w.imageFlags, s, !1, 0, 0), "THREE.BasisTextureLoader: transcodeImage() failed for level " + x.index + "."), A.push({
												data: S,
												width: x.width,
												height: x.height
											})
										}
									} finally {
										g.delete()
									}
								} else
									for (b = 0; b < e.levels.length; b++) {
										x = e.levels[b], C = p(o, x.width, x.height), S = new Uint8Array(C);
										h(i.transcodeUASTCImage(o, S, C / m, x.data, u(o, x.width), d(o, x.height), x.width, x.height, x.index, 0, x.data.byteLength, 0, s, !1, 0, 0, -1, -1), "THREE.BasisTextureLoader: transcodeUASTCImage() failed for level " + x.index + "."), A.push({
											data: S,
											width: x.width,
											height: x.height
										})
									}
								return {
									width: n,
									height: r,
									hasAlpha: s,
									mipmaps: A,
									format: l
								}
							}(s.taskConfig) : function(e) {
								var t = new i.BasisFile(new Uint8Array(e)),
									n = t.isUASTC() ? a.UASTC_4x4 : a.ETC1S,
									r = t.getImageWidth(0, 0),
									s = t.getImageHeight(0, 0),
									o = t.getNumLevels(0),
									l = t.getHasAlpha();

								function h() {
									t.close(), t.delete()
								}
								var {
									transcoderFormat: u,
									engineFormat: d
								} = c(n, r, s, l);
								if (!r || !s || !o) throw h(), new Error("THREE.BasisTextureLoader:\tInvalid texture");
								if (!t.startTranscoding()) throw h(), new Error("THREE.BasisTextureLoader: .startTranscoding failed");
								for (var p = [], m = 0; m < o; m++) {
									var A = t.getImageWidth(0, m),
										g = t.getImageHeight(0, m),
										f = new Uint8Array(t.getImageTranscodedSizeInBytes(0, m, u));
									if (!t.transcodeImage(f, 0, m, u, 0, l)) throw h(), new Error("THREE.BasisTextureLoader: .transcodeImage failed.");
									p.push({
										data: f,
										width: A,
										height: g
									})
								}
								return h(), {
									width: r,
									height: s,
									hasAlpha: l,
									mipmaps: p,
									format: d
								}
							}(s.buffers[0]), l = [], m = 0; m < r.length; ++m) l.push(r[m].data.buffer);
						self.postMessage({
							type: "transcode",
							id: s.id,
							width: e,
							height: t,
							hasAlpha: n,
							mipmaps: r,
							format: o
						}, l)
					} catch (e) {
						console.error(e), self.postMessage({
							type: "error",
							id: s.id,
							error: e.message
						})
					}
				})
		}
	};
	var s = [{
			if: "astcSupported",
			basisFormat: [a.UASTC_4x4],
			transcoderFormat: [r.ASTC_4x4, r.ASTC_4x4],
			engineFormat: [n.RGBA_ASTC_4x4_Format, n.RGBA_ASTC_4x4_Format],
			priorityETC1S: 1 / 0,
			priorityUASTC: 1,
			needsPowerOfTwo: !1
		}, {
			if: "bptcSupported",
			basisFormat: [a.ETC1S, a.UASTC_4x4],
			transcoderFormat: [r.BC7_M5, r.BC7_M5],
			engineFormat: [n.RGBA_BPTC_Format, n.RGBA_BPTC_Format],
			priorityETC1S: 3,
			priorityUASTC: 2,
			needsPowerOfTwo: !1
		}, {
			if: "dxtSupported",
			basisFormat: [a.ETC1S, a.UASTC_4x4],
			transcoderFormat: [r.BC1, r.BC3],
			engineFormat: [n.RGB_S3TC_DXT1_Format, n.RGBA_S3TC_DXT5_Format],
			priorityETC1S: 4,
			priorityUASTC: 5,
			needsPowerOfTwo: !1
		}, {
			if: "etc2Supported",
			basisFormat: [a.ETC1S, a.UASTC_4x4],
			transcoderFormat: [r.ETC1, r.ETC2],
			engineFormat: [n.RGB_ETC2_Format, n.RGBA_ETC2_EAC_Format],
			priorityETC1S: 1,
			priorityUASTC: 3,
			needsPowerOfTwo: !1
		}, {
			if: "etc1Supported",
			basisFormat: [a.ETC1S, a.UASTC_4x4],
			transcoderFormat: [r.ETC1, r.ETC1],
			engineFormat: [n.RGB_ETC1_Format, n.RGB_ETC1_Format],
			priorityETC1S: 2,
			priorityUASTC: 4,
			needsPowerOfTwo: !1
		}, {
			if: "pvrtcSupported",
			basisFormat: [a.ETC1S, a.UASTC_4x4],
			transcoderFormat: [r.PVRTC1_4_RGB, r.PVRTC1_4_RGBA],
			engineFormat: [n.RGB_PVRTC_4BPPV1_Format, n.RGBA_PVRTC_4BPPV1_Format],
			priorityETC1S: 5,
			priorityUASTC: 6,
			needsPowerOfTwo: !0
		}],
		o = s.sort((function(e, t) {
			return e.priorityETC1S - t.priorityETC1S
		})),
		l = s.sort((function(e, t) {
			return e.priorityUASTC - t.priorityUASTC
		}));

	function c(t, i, s, c) {
		for (var h = t === a.ETC1S ? o : l, u = 0; u < h.length; u++) {
			var d = h[u];
			if (e[d.if] && (d.basisFormat.includes(t) && (!d.needsPowerOfTwo || m(i) && m(s)))) return {
				transcoderFormat: d.transcoderFormat[c ? 1 : 0],
				engineFormat: d.engineFormat[c ? 1 : 0]
			}
		}
		return console.warn("THREE.BasisTextureLoader: No suitable compressed texture format found. Decoding to RGBA32."), {
			transcoderFormat: r.RGBA32,
			engineFormat: n.RGBAFormat
		}
	}

	function h(e, t) {
		if (!e) throw new Error(t)
	}

	function u(e, t) {
		return Math.ceil(t / i.getFormatBlockWidth(e))
	}

	function d(e, t) {
		return Math.ceil(t / i.getFormatBlockHeight(e))
	}

	function p(e, t, n) {
		var a = i.getBytesPerBlockOrPixel(e);
		if (i.formatIsUncompressed(e)) return t * n * a;
		if (e === r.PVRTC1_4_RGB || e === r.PVRTC1_4_RGBA) {
			var s = t + 3 & -4,
				o = n + 3 & -4;
			return (Math.max(8, s) * Math.max(8, o) * 4 + 7) / 8
		}
		return u(e, t) * d(e, n) * a
	}

	function m(e) {
		return e <= 2 || 0 == (e & e - 1) && 0 !== e
	}
};
const importObject = {
	env: {
		emscripten_notify_memory_growth: function(e) {
			heap = new Uint8Array(instance.exports.memory.buffer)
		}
	}
};
class ZSTDDecoder {
	init() {
		return init || (init = fetch("data:application/wasm;base64," + wasm).then(e => e.arrayBuffer()).then(e => WebAssembly.instantiate(e, importObject)).then(e => {
			instance = e.instance, importObject.env.emscripten_notify_memory_growth(0)
		})), init
	}
	decode(e, t = 0) {
		const i = e.byteLength,
			n = instance.exports.malloc(i);
		heap.set(e, n), t = t || Number(instance.exports.ZSTD_findDecompressedSize(n, i));
		const r = instance.exports.malloc(t),
			a = instance.exports.ZSTD_decompress(r, t, n, i),
			s = heap.slice(r, r + a);
		return instance.exports.free(n), instance.exports.free(r), s
	}
}
const wasm = "AGFzbQEAAAABpQEVYAF/AX9gAn9/AGADf39/AX9gBX9/f39/AX9gAX8AYAJ/fwF/YAR/f39/AX9gA39/fwBgBn9/f39/fwF/YAd/f39/f39/AX9gAn9/AX5gAn5+AX5gAABgBX9/f39/AGAGf39/f39/AGAIf39/f39/f38AYAl/f39/f39/f38AYAABf2AIf39/f39/f38Bf2ANf39/f39/f39/f39/fwF/YAF/AX4CJwEDZW52H2Vtc2NyaXB0ZW5fbm90aWZ5X21lbW9yeV9ncm93dGgABANpaAEFAAAFAgEFCwACAQABAgIFBQcAAwABDgsBAQcAEhMHAAUBDAQEAAANBwQCAgYCBAgDAwMDBgEACQkHBgICAAYGAgQUBwYGAwIGAAMCAQgBBwUGCgoEEQAEBAEIAwgDBQgDEA8IAAcABAUBcAECAgUEAQCAAgYJAX8BQaCgwAILB2AHBm1lbW9yeQIABm1hbGxvYwAoBGZyZWUAJgxaU1REX2lzRXJyb3IAaBlaU1REX2ZpbmREZWNvbXByZXNzZWRTaXplAFQPWlNURF9kZWNvbXByZXNzAEoGX3N0YXJ0ACQJBwEAQQELASQKussBaA8AIAAgACgCBCABajYCBAsZACAAKAIAIAAoAgRBH3F0QQAgAWtBH3F2CwgAIABBiH9LC34BBH9BAyEBIAAoAgQiA0EgTQRAIAAoAggiASAAKAIQTwRAIAAQDQ8LIAAoAgwiAiABRgRAQQFBAiADQSBJGw8LIAAgASABIAJrIANBA3YiBCABIARrIAJJIgEbIgJrIgQ2AgggACADIAJBA3RrNgIEIAAgBCgAADYCAAsgAQsUAQF/IAAgARACIQIgACABEAEgAgv3AQECfyACRQRAIABCADcCACAAQQA2AhAgAEIANwIIQbh/DwsgACABNgIMIAAgAUEEajYCECACQQRPBEAgACABIAJqIgFBfGoiAzYCCCAAIAMoAAA2AgAgAUF/ai0AACIBBEAgAEEIIAEQFGs2AgQgAg8LIABBADYCBEF/DwsgACABNgIIIAAgAS0AACIDNgIAIAJBfmoiBEEBTQRAIARBAWtFBEAgACABLQACQRB0IANyIgM2AgALIAAgAS0AAUEIdCADajYCAAsgASACakF/ai0AACIBRQRAIABBADYCBEFsDwsgAEEoIAEQFCACQQN0ams2AgQgAgsWACAAIAEpAAA3AAAgACABKQAINwAICy8BAX8gAUECdEGgHWooAgAgACgCAEEgIAEgACgCBGprQR9xdnEhAiAAIAEQASACCyEAIAFCz9bTvtLHq9lCfiAAfEIfiUKHla+vmLbem55/fgsdAQF/IAAoAgggACgCDEYEfyAAKAIEQSBGBUEACwuCBAEDfyACQYDAAE8EQCAAIAEgAhBnIAAPCyAAIAJqIQMCQCAAIAFzQQNxRQRAAkAgAkEBSARAIAAhAgwBCyAAQQNxRQRAIAAhAgwBCyAAIQIDQCACIAEtAAA6AAAgAUEBaiEBIAJBAWoiAiADTw0BIAJBA3ENAAsLAkAgA0F8cSIEQcAASQ0AIAIgBEFAaiIFSw0AA0AgAiABKAIANgIAIAIgASgCBDYCBCACIAEoAgg2AgggAiABKAIMNgIMIAIgASgCEDYCECACIAEoAhQ2AhQgAiABKAIYNgIYIAIgASgCHDYCHCACIAEoAiA2AiAgAiABKAIkNgIkIAIgASgCKDYCKCACIAEoAiw2AiwgAiABKAIwNgIwIAIgASgCNDYCNCACIAEoAjg2AjggAiABKAI8NgI8IAFBQGshASACQUBrIgIgBU0NAAsLIAIgBE8NAQNAIAIgASgCADYCACABQQRqIQEgAkEEaiICIARJDQALDAELIANBBEkEQCAAIQIMAQsgA0F8aiIEIABJBEAgACECDAELIAAhAgNAIAIgAS0AADoAACACIAEtAAE6AAEgAiABLQACOgACIAIgAS0AAzoAAyABQQRqIQEgAkEEaiICIARNDQALCyACIANJBEADQCACIAEtAAA6AAAgAUEBaiEBIAJBAWoiAiADRw0ACwsgAAsMACAAIAEpAAA3AAALQQECfyAAKAIIIgEgACgCEEkEQEEDDwsgACAAKAIEIgJBB3E2AgQgACABIAJBA3ZrIgE2AgggACABKAAANgIAQQALDAAgACABKAIANgAAC/cCAQJ/AkAgACABRg0AAkAgASACaiAASwRAIAAgAmoiBCABSw0BCyAAIAEgAhALDwsgACABc0EDcSEDAkACQCAAIAFJBEAgAwRAIAAhAwwDCyAAQQNxRQRAIAAhAwwCCyAAIQMDQCACRQ0EIAMgAS0AADoAACABQQFqIQEgAkF/aiECIANBAWoiA0EDcQ0ACwwBCwJAIAMNACAEQQNxBEADQCACRQ0FIAAgAkF/aiICaiIDIAEgAmotAAA6AAAgA0EDcQ0ACwsgAkEDTQ0AA0AgACACQXxqIgJqIAEgAmooAgA2AgAgAkEDSw0ACwsgAkUNAgNAIAAgAkF/aiICaiABIAJqLQAAOgAAIAINAAsMAgsgAkEDTQ0AIAIhBANAIAMgASgCADYCACABQQRqIQEgA0EEaiEDIARBfGoiBEEDSw0ACyACQQNxIQILIAJFDQADQCADIAEtAAA6AAAgA0EBaiEDIAFBAWohASACQX9qIgINAAsLIAAL8wICAn8BfgJAIAJFDQAgACACaiIDQX9qIAE6AAAgACABOgAAIAJBA0kNACADQX5qIAE6AAAgACABOgABIANBfWogAToAACAAIAE6AAIgAkEHSQ0AIANBfGogAToAACAAIAE6AAMgAkEJSQ0AIABBACAAa0EDcSIEaiIDIAFB/wFxQYGChAhsIgE2AgAgAyACIARrQXxxIgRqIgJBfGogATYCACAEQQlJDQAgAyABNgIIIAMgATYCBCACQXhqIAE2AgAgAkF0aiABNgIAIARBGUkNACADIAE2AhggAyABNgIUIAMgATYCECADIAE2AgwgAkFwaiABNgIAIAJBbGogATYCACACQWhqIAE2AgAgAkFkaiABNgIAIAQgA0EEcUEYciIEayICQSBJDQAgAa0iBUIghiAFhCEFIAMgBGohAQNAIAEgBTcDGCABIAU3AxAgASAFNwMIIAEgBTcDACABQSBqIQEgAkFgaiICQR9LDQALCyAACy8BAn8gACgCBCAAKAIAQQJ0aiICLQACIQMgACACLwEAIAEgAi0AAxAIajYCACADCy8BAn8gACgCBCAAKAIAQQJ0aiICLQACIQMgACACLwEAIAEgAi0AAxAFajYCACADCx8AIAAgASACKAIEEAg2AgAgARAEGiAAIAJBCGo2AgQLCAAgAGdBH3MLugUBDX8jAEEQayIKJAACfyAEQQNNBEAgCkEANgIMIApBDGogAyAEEAsaIAAgASACIApBDGpBBBAVIgBBbCAAEAMbIAAgACAESxsMAQsgAEEAIAEoAgBBAXRBAmoQECENQVQgAygAACIGQQ9xIgBBCksNABogAiAAQQVqNgIAIAMgBGoiAkF8aiEMIAJBeWohDiACQXtqIRAgAEEGaiELQQQhBSAGQQR2IQRBICAAdCIAQQFyIQkgASgCACEPQQAhAiADIQYCQANAIAlBAkggAiAPS3JFBEAgAiEHAkAgCARAA0AgBEH//wNxQf//A0YEQCAHQRhqIQcgBiAQSQR/IAZBAmoiBigAACAFdgUgBUEQaiEFIARBEHYLIQQMAQsLA0AgBEEDcSIIQQNGBEAgBUECaiEFIARBAnYhBCAHQQNqIQcMAQsLIAcgCGoiByAPSw0EIAVBAmohBQNAIAIgB0kEQCANIAJBAXRqQQA7AQAgAkEBaiECDAELCyAGIA5LQQAgBiAFQQN1aiIHIAxLG0UEQCAHKAAAIAVBB3EiBXYhBAwCCyAEQQJ2IQQLIAYhBwsCfyALQX9qIAQgAEF/anEiBiAAQQF0QX9qIgggCWsiEUkNABogBCAIcSIEQQAgESAEIABIG2shBiALCyEIIA0gAkEBdGogBkF/aiIEOwEAIAlBASAGayAEIAZBAUgbayEJA0AgCSAASARAIABBAXUhACALQX9qIQsMAQsLAn8gByAOS0EAIAcgBSAIaiIFQQN1aiIGIAxLG0UEQCAFQQdxDAELIAUgDCIGIAdrQQN0awshBSACQQFqIQIgBEUhCCAGKAAAIAVBH3F2IQQMAQsLQWwgCUEBRyAFQSBKcg0BGiABIAJBf2o2AgAgBiAFQQdqQQN1aiADawwBC0FQCyEAIApBEGokACAACwkAQQFBBSAAGwsMACAAIAEoAAA2AAALqgMBCn8jAEHwAGsiCiQAIAJBAWohDiAAQQhqIQtBgIAEIAVBf2p0QRB1IQxBACECQQEhBkEBIAV0IglBf2oiDyEIA0AgAiAORkUEQAJAIAEgAkEBdCINai8BACIHQf//A0YEQCALIAhBA3RqIAI2AgQgCEF/aiEIQQEhBwwBCyAGQQAgDCAHQRB0QRB1ShshBgsgCiANaiAHOwEAIAJBAWohAgwBCwsgACAFNgIEIAAgBjYCACAJQQN2IAlBAXZqQQNqIQxBACEAQQAhBkEAIQIDQCAGIA5GBEADQAJAIAAgCUYNACAKIAsgAEEDdGoiASgCBCIGQQF0aiICIAIvAQAiAkEBajsBACABIAUgAhAUayIIOgADIAEgAiAIQf8BcXQgCWs7AQAgASAEIAZBAnQiAmooAgA6AAIgASACIANqKAIANgIEIABBAWohAAwBCwsFIAEgBkEBdGouAQAhDUEAIQcDQCAHIA1ORQRAIAsgAkEDdGogBjYCBANAIAIgDGogD3EiAiAISw0ACyAHQQFqIQcMAQsLIAZBAWohBgwBCwsgCkHwAGokAAsjAEIAIAEQCSAAhUKHla+vmLbem55/fkLj3MqV/M7y9YV/fAsQACAAQn43AwggACABNgIACyQBAX8gAARAIAEoAgQiAgRAIAEoAgggACACEQEADwsgABAmCwsfACAAIAEgAi8BABAINgIAIAEQBBogACACQQRqNgIEC0oBAX9BoCAoAgAiASAAaiIAQX9MBEBBiCBBMDYCAEF/DwsCQCAAPwBBEHRNDQAgABBmDQBBiCBBMDYCAEF/DwtBoCAgADYCACABC9cBAQh/Qbp/IQoCQCACKAIEIgggAigCACIJaiIOIAEgAGtLDQBBbCEKIAkgBCADKAIAIgtrSw0AIAAgCWoiBCACKAIIIgxrIQ0gACABQWBqIg8gCyAJQQAQKSADIAkgC2o2AgACQAJAIAwgBCAFa00EQCANIQUMAQsgDCAEIAZrSw0CIAcgDSAFayIAaiIBIAhqIAdNBEAgBCABIAgQDxoMAgsgBCABQQAgAGsQDyEBIAIgACAIaiIINgIEIAEgAGshBAsgBCAPIAUgCEEBECkLIA4hCgsgCgubAgEBfyMAQYABayINJAAgDSADNgJ8AkAgAkEDSwRAQX8hCQwBCwJAAkACQAJAIAJBAWsOAwADAgELIAZFBEBBuH8hCQwEC0FsIQkgBS0AACICIANLDQMgACAHIAJBAnQiAmooAgAgAiAIaigCABA7IAEgADYCAEEBIQkMAwsgASAJNgIAQQAhCQwCCyAKRQRAQWwhCQwCC0EAIQkgC0UgDEEZSHINAUEIIAR0QQhqIQBBACECA0AgAiAATw0CIAJBQGshAgwAAAsAC0FsIQkgDSANQfwAaiANQfgAaiAFIAYQFSICEAMNACANKAJ4IgMgBEsNACAAIA0gDSgCfCAHIAggAxAYIAEgADYCACACIQkLIA1BgAFqJAAgCQsLACAAIAEgAhALGgsQACAALwAAIAAtAAJBEHRyCy8AAn9BuH8gAUEISQ0AGkFyIAAoAAQiAEF3Sw0AGkG4fyAAQQhqIgAgACABSxsLCwkAIAAgATsAAAsDAAELigYBBX8gACAAKAIAIgVBfnE2AgBBACAAIAVBAXZqQYQgKAIAIgQgAEYbIQECQAJAIAAoAgQiAkUNACACKAIAIgNBAXENACACQQhqIgUgA0EBdkF4aiIDQQggA0EISxtnQR9zQQJ0QYAfaiIDKAIARgRAIAMgAigCDDYCAAsgAigCCCIDBEAgAyACKAIMNgIECyACKAIMIgMEQCADIAIoAgg2AgALIAIgAigCACAAKAIAQX5xajYCAEGEICEAAkACQCABRQ0AIAEgAjYCBCABKAIAIgNBAXENASADQQF2QXhqIgNBCCADQQhLG2dBH3NBAnRBgB9qIgMoAgAgAUEIakYEQCADIAEoAgw2AgALIAEoAggiAwRAIAMgASgCDDYCBAsgASgCDCIDBEAgAyABKAIINgIAQYQgKAIAIQQLIAIgAigCACABKAIAQX5xajYCACABIARGDQAgASABKAIAQQF2akEEaiEACyAAIAI2AgALIAIoAgBBAXZBeGoiAEEIIABBCEsbZ0Efc0ECdEGAH2oiASgCACEAIAEgBTYCACACIAA2AgwgAkEANgIIIABFDQEgACAFNgIADwsCQCABRQ0AIAEoAgAiAkEBcQ0AIAJBAXZBeGoiAkEIIAJBCEsbZ0Efc0ECdEGAH2oiAigCACABQQhqRgRAIAIgASgCDDYCAAsgASgCCCICBEAgAiABKAIMNgIECyABKAIMIgIEQCACIAEoAgg2AgBBhCAoAgAhBAsgACAAKAIAIAEoAgBBfnFqIgI2AgACQCABIARHBEAgASABKAIAQQF2aiAANgIEIAAoAgAhAgwBC0GEICAANgIACyACQQF2QXhqIgFBCCABQQhLG2dBH3NBAnRBgB9qIgIoAgAhASACIABBCGoiAjYCACAAIAE2AgwgAEEANgIIIAFFDQEgASACNgIADwsgBUEBdkF4aiIBQQggAUEISxtnQR9zQQJ0QYAfaiICKAIAIQEgAiAAQQhqIgI2AgAgACABNgIMIABBADYCCCABRQ0AIAEgAjYCAAsLDgAgAARAIABBeGoQJQsLgAIBA38CQCAAQQ9qQXhxQYQgKAIAKAIAQQF2ayICEB1Bf0YNAAJAQYQgKAIAIgAoAgAiAUEBcQ0AIAFBAXZBeGoiAUEIIAFBCEsbZ0Efc0ECdEGAH2oiASgCACAAQQhqRgRAIAEgACgCDDYCAAsgACgCCCIBBEAgASAAKAIMNgIECyAAKAIMIgFFDQAgASAAKAIINgIAC0EBIQEgACAAKAIAIAJBAXRqIgI2AgAgAkEBcQ0AIAJBAXZBeGoiAkEIIAJBCEsbZ0Efc0ECdEGAH2oiAygCACECIAMgAEEIaiIDNgIAIAAgAjYCDCAAQQA2AgggAkUNACACIAM2AgALIAELtwIBA38CQAJAIABBASAAGyICEDgiAA0AAkACQEGEICgCACIARQ0AIAAoAgAiA0EBcQ0AIAAgA0EBcjYCACADQQF2QXhqIgFBCCABQQhLG2dBH3NBAnRBgB9qIgEoAgAgAEEIakYEQCABIAAoAgw2AgALIAAoAggiAQRAIAEgACgCDDYCBAsgACgCDCIBBEAgASAAKAIINgIACyACECchAkEAIQFBhCAoAgAhACACDQEgACAAKAIAQX5xNgIAQQAPCyACQQ9qQXhxIgMQHSICQX9GDQIgAkEHakF4cSIAIAJHBEAgACACaxAdQX9GDQMLAkBBhCAoAgAiAUUEQEGAICAANgIADAELIAAgATYCBAtBhCAgADYCACAAIANBAXRBAXI2AgAMAQsgAEUNAQsgAEEIaiEBCyABC7kDAQJ/IAAgA2ohBQJAIANBB0wEQANAIAAgBU8NAiAAIAItAAA6AAAgAEEBaiEAIAJBAWohAgwAAAsACyAEQQFGBEACQCAAIAJrIgZBB00EQCAAIAItAAA6AAAgACACLQABOgABIAAgAi0AAjoAAiAAIAItAAM6AAMgAEEEaiACIAZBAnQiBkHAHmooAgBqIgIQFyACIAZB4B5qKAIAayECDAELIAAgAhAMCyACQQhqIQIgAEEIaiEACwJAAkACQAJAIAUgAU0EQCAAIANqIQEgBEEBRyAAIAJrQQ9Kcg0BA0AgACACEAwgAkEIaiECIABBCGoiACABSQ0ACwwFCyAAIAFLBEAgACEBDAQLIARBAUcgACACa0EPSnINASAAIQMgAiEEA0AgAyAEEAwgBEEIaiEEIANBCGoiAyABSQ0ACwwCCwNAIAAgAhAHIAJBEGohAiAAQRBqIgAgAUkNAAsMAwsgACEDIAIhBANAIAMgBBAHIARBEGohBCADQRBqIgMgAUkNAAsLIAIgASAAa2ohAgsDQCABIAVPDQEgASACLQAAOgAAIAFBAWohASACQQFqIQIMAAALAAsLQQECfyAAIAAoArjgASIDNgLE4AEgACgCvOABIQQgACABNgK84AEgACABIAJqNgK44AEgACABIAQgA2tqNgLA4AELpgEBAX8gACAAKALs4QEQFjYCyOABIABCADcD+OABIABCADcDuOABIABBwOABakIANwMAIABBqNAAaiIBQYyAgOAANgIAIABBADYCmOIBIABCADcDiOEBIABCAzcDgOEBIABBrNABakHgEikCADcCACAAQbTQAWpB6BIoAgA2AgAgACABNgIMIAAgAEGYIGo2AgggACAAQaAwajYCBCAAIABBEGo2AgALYQEBf0G4fyEDAkAgAUEDSQ0AIAIgABAhIgFBA3YiADYCCCACIAFBAXE2AgQgAiABQQF2QQNxIgM2AgACQCADQX9qIgFBAksNAAJAIAFBAWsOAgEAAgtBbA8LIAAhAwsgAwsMACAAIAEgAkEAEC4LiAQCA38CfiADEBYhBCAAQQBBKBAQIQAgBCACSwRAIAQPCyABRQRAQX8PCwJAAkAgA0EBRg0AIAEoAAAiBkGo6r5pRg0AQXYhAyAGQXBxQdDUtMIBRw0BQQghAyACQQhJDQEgAEEAQSgQECEAIAEoAAQhASAAQQE2AhQgACABrTcDAEEADwsgASACIAMQLyIDIAJLDQAgACADNgIYQXIhAyABIARqIgVBf2otAAAiAkEIcQ0AIAJBIHEiBkUEQEFwIQMgBS0AACIFQacBSw0BIAVBB3GtQgEgBUEDdkEKaq2GIgdCA4h+IAd8IQggBEEBaiEECyACQQZ2IQMgAkECdiEFAkAgAkEDcUF/aiICQQJLBEBBACECDAELAkACQAJAIAJBAWsOAgECAAsgASAEai0AACECIARBAWohBAwCCyABIARqLwAAIQIgBEECaiEEDAELIAEgBGooAAAhAiAEQQRqIQQLIAVBAXEhBQJ+AkACQAJAIANBf2oiA0ECTQRAIANBAWsOAgIDAQtCfyAGRQ0DGiABIARqMQAADAMLIAEgBGovAACtQoACfAwCCyABIARqKAAArQwBCyABIARqKQAACyEHIAAgBTYCICAAIAI2AhwgACAHNwMAQQAhAyAAQQA2AhQgACAHIAggBhsiBzcDCCAAIAdCgIAIIAdCgIAIVBs+AhALIAMLWwEBf0G4fyEDIAIQFiICIAFNBH8gACACakF/ai0AACIAQQNxQQJ0QaAeaigCACACaiAAQQZ2IgFBAnRBsB5qKAIAaiAAQSBxIgBFaiABRSAAQQV2cWoFQbh/CwsdACAAKAKQ4gEQWiAAQQA2AqDiASAAQgA3A5DiAQu1AwEFfyMAQZACayIKJABBuH8hBgJAIAVFDQAgBCwAACIIQf8BcSEHAkAgCEF/TARAIAdBgn9qQQF2IgggBU8NAkFsIQYgB0GBf2oiBUGAAk8NAiAEQQFqIQdBACEGA0AgBiAFTwRAIAUhBiAIIQcMAwUgACAGaiAHIAZBAXZqIgQtAABBBHY6AAAgACAGQQFyaiAELQAAQQ9xOgAAIAZBAmohBgwBCwAACwALIAcgBU8NASAAIARBAWogByAKEFMiBhADDQELIAYhBEEAIQYgAUEAQTQQECEJQQAhBQNAIAQgBkcEQCAAIAZqIggtAAAiAUELSwRAQWwhBgwDBSAJIAFBAnRqIgEgASgCAEEBajYCACAGQQFqIQZBASAILQAAdEEBdSAFaiEFDAILAAsLQWwhBiAFRQ0AIAUQFEEBaiIBQQxLDQAgAyABNgIAQQFBASABdCAFayIDEBQiAXQgA0cNACAAIARqIAFBAWoiADoAACAJIABBAnRqIgAgACgCAEEBajYCACAJKAIEIgBBAkkgAEEBcXINACACIARBAWo2AgAgB0EBaiEGCyAKQZACaiQAIAYLxhEBDH8jAEHwAGsiBSQAQWwhCwJAIANBCkkNACACLwAAIQogAi8AAiEJIAIvAAQhByAFQQhqIAQQDgJAIAMgByAJIApqakEGaiIMSQ0AIAUtAAohCCAFQdgAaiACQQZqIgIgChAGIgsQAw0BIAVBQGsgAiAKaiICIAkQBiILEAMNASAFQShqIAIgCWoiAiAHEAYiCxADDQEgBUEQaiACIAdqIAMgDGsQBiILEAMNASAAIAFqIg9BfWohECAEQQRqIQZBASELIAAgAUEDakECdiIDaiIMIANqIgIgA2oiDiEDIAIhBCAMIQcDQCALIAMgEElxBEAgACAGIAVB2ABqIAgQAkECdGoiCS8BADsAACAFQdgAaiAJLQACEAEgCS0AAyELIAcgBiAFQUBrIAgQAkECdGoiCS8BADsAACAFQUBrIAktAAIQASAJLQADIQogBCAGIAVBKGogCBACQQJ0aiIJLwEAOwAAIAVBKGogCS0AAhABIAktAAMhCSADIAYgBUEQaiAIEAJBAnRqIg0vAQA7AAAgBUEQaiANLQACEAEgDS0AAyENIAAgC2oiCyAGIAVB2ABqIAgQAkECdGoiAC8BADsAACAFQdgAaiAALQACEAEgAC0AAyEAIAcgCmoiCiAGIAVBQGsgCBACQQJ0aiIHLwEAOwAAIAVBQGsgBy0AAhABIActAAMhByAEIAlqIgkgBiAFQShqIAgQAkECdGoiBC8BADsAACAFQShqIAQtAAIQASAELQADIQQgAyANaiIDIAYgBUEQaiAIEAJBAnRqIg0vAQA7AAAgBUEQaiANLQACEAEgACALaiEAIAcgCmohByAEIAlqIQQgAyANLQADaiEDIAVB2ABqEA0gBUFAaxANciAFQShqEA1yIAVBEGoQDXJFIQsMAQsLIAQgDksgByACS3INAEFsIQsgACAMSw0BIAxBfWohCQNAQQAgACAJSSAFQdgAahAEGwRAIAAgBiAFQdgAaiAIEAJBAnRqIgovAQA7AAAgBUHYAGogCi0AAhABIAAgCi0AA2oiACAGIAVB2ABqIAgQAkECdGoiCi8BADsAACAFQdgAaiAKLQACEAEgACAKLQADaiEADAEFIAxBfmohCgNAIAVB2ABqEAQgACAKS3JFBEAgACAGIAVB2ABqIAgQAkECdGoiCS8BADsAACAFQdgAaiAJLQACEAEgACAJLQADaiEADAELCwNAIAAgCk0EQCAAIAYgBUHYAGogCBACQQJ0aiIJLwEAOwAAIAVB2ABqIAktAAIQASAAIAktAANqIQAMAQsLAkAgACAMTw0AIAAgBiAFQdgAaiAIEAIiAEECdGoiDC0AADoAACAMLQADQQFGBEAgBUHYAGogDC0AAhABDAELIAUoAlxBH0sNACAFQdgAaiAGIABBAnRqLQACEAEgBSgCXEEhSQ0AIAVBIDYCXAsgAkF9aiEMA0BBACAHIAxJIAVBQGsQBBsEQCAHIAYgBUFAayAIEAJBAnRqIgAvAQA7AAAgBUFAayAALQACEAEgByAALQADaiIAIAYgBUFAayAIEAJBAnRqIgcvAQA7AAAgBUFAayAHLQACEAEgACAHLQADaiEHDAEFIAJBfmohDANAIAVBQGsQBCAHIAxLckUEQCAHIAYgBUFAayAIEAJBAnRqIgAvAQA7AAAgBUFAayAALQACEAEgByAALQADaiEHDAELCwNAIAcgDE0EQCAHIAYgBUFAayAIEAJBAnRqIgAvAQA7AAAgBUFAayAALQACEAEgByAALQADaiEHDAELCwJAIAcgAk8NACAHIAYgBUFAayAIEAIiAEECdGoiAi0AADoAACACLQADQQFGBEAgBUFAayACLQACEAEMAQsgBSgCREEfSw0AIAVBQGsgBiAAQQJ0ai0AAhABIAUoAkRBIUkNACAFQSA2AkQLIA5BfWohAgNAQQAgBCACSSAFQShqEAQbBEAgBCAGIAVBKGogCBACQQJ0aiIALwEAOwAAIAVBKGogAC0AAhABIAQgAC0AA2oiACAGIAVBKGogCBACQQJ0aiIELwEAOwAAIAVBKGogBC0AAhABIAAgBC0AA2ohBAwBBSAOQX5qIQIDQCAFQShqEAQgBCACS3JFBEAgBCAGIAVBKGogCBACQQJ0aiIALwEAOwAAIAVBKGogAC0AAhABIAQgAC0AA2ohBAwBCwsDQCAEIAJNBEAgBCAGIAVBKGogCBACQQJ0aiIALwEAOwAAIAVBKGogAC0AAhABIAQgAC0AA2ohBAwBCwsCQCAEIA5PDQAgBCAGIAVBKGogCBACIgBBAnRqIgItAAA6AAAgAi0AA0EBRgRAIAVBKGogAi0AAhABDAELIAUoAixBH0sNACAFQShqIAYgAEECdGotAAIQASAFKAIsQSFJDQAgBUEgNgIsCwNAQQAgAyAQSSAFQRBqEAQbBEAgAyAGIAVBEGogCBACQQJ0aiIALwEAOwAAIAVBEGogAC0AAhABIAMgAC0AA2oiACAGIAVBEGogCBACQQJ0aiICLwEAOwAAIAVBEGogAi0AAhABIAAgAi0AA2ohAwwBBSAPQX5qIQIDQCAFQRBqEAQgAyACS3JFBEAgAyAGIAVBEGogCBACQQJ0aiIALwEAOwAAIAVBEGogAC0AAhABIAMgAC0AA2ohAwwBCwsDQCADIAJNBEAgAyAGIAVBEGogCBACQQJ0aiIALwEAOwAAIAVBEGogAC0AAhABIAMgAC0AA2ohAwwBCwsCQCADIA9PDQAgAyAGIAVBEGogCBACIgBBAnRqIgItAAA6AAAgAi0AA0EBRgRAIAVBEGogAi0AAhABDAELIAUoAhRBH0sNACAFQRBqIAYgAEECdGotAAIQASAFKAIUQSFJDQAgBUEgNgIUCyABQWwgBUHYAGoQCiAFQUBrEApxIAVBKGoQCnEgBUEQahAKcRshCwwJCwAACwALAAALAAsAAAsACwAACwALQWwhCwsgBUHwAGokACALC7UEAQ5/IwBBEGsiBiQAIAZBBGogABAOQVQhBQJAIARB3AtJDQAgBi0ABCEHIANB8ARqQQBB7AAQECEIIAdBDEsNACADQdwJaiIJIAggBkEIaiAGQQxqIAEgAhAxIhAQA0UEQCAGKAIMIgQgB0sNASADQdwFaiEPIANBpAVqIREgAEEEaiESIANBqAVqIQEgBCEFA0AgBSICQX9qIQUgCCACQQJ0aigCAEUNAAsgAkEBaiEOQQEhBQNAIAUgDk9FBEAgCCAFQQJ0IgtqKAIAIQwgASALaiAKNgIAIAVBAWohBSAKIAxqIQoMAQsLIAEgCjYCAEEAIQUgBigCCCELA0AgBSALRkUEQCABIAUgCWotAAAiDEECdGoiDSANKAIAIg1BAWo2AgAgDyANQQF0aiINIAw6AAEgDSAFOgAAIAVBAWohBQwBCwtBACEBIANBADYCqAUgBEF/cyAHaiEJQQEhBQNAIAUgDk9FBEAgCCAFQQJ0IgtqKAIAIQwgAyALaiABNgIAIAwgBSAJanQgAWohASAFQQFqIQUMAQsLIAcgBEEBaiIBIAJrIgRrQQFqIQgDQEEBIQUgBCAIT0UEQANAIAUgDk9FBEAgBUECdCIJIAMgBEE0bGpqIAMgCWooAgAgBHY2AgAgBUEBaiEFDAELCyAEQQFqIQQMAQsLIBIgByAPIAogESADIAIgARBkIAZBAToABSAGIAc6AAYgACAGKAIENgIACyAQIQULIAZBEGokACAFC8ENAQt/IwBB8ABrIgUkAEFsIQkCQCADQQpJDQAgAi8AACEKIAIvAAIhDCACLwAEIQYgBUEIaiAEEA4CQCADIAYgCiAMampBBmoiDUkNACAFLQAKIQcgBUHYAGogAkEGaiICIAoQBiIJEAMNASAFQUBrIAIgCmoiAiAMEAYiCRADDQEgBUEoaiACIAxqIgIgBhAGIgkQAw0BIAVBEGogAiAGaiADIA1rEAYiCRADDQEgACABaiIOQX1qIQ8gBEEEaiEGQQEhCSAAIAFBA2pBAnYiAmoiCiACaiIMIAJqIg0hAyAMIQQgCiECA0AgCSADIA9JcQRAIAYgBUHYAGogBxACQQF0aiIILQAAIQsgBUHYAGogCC0AARABIAAgCzoAACAGIAVBQGsgBxACQQF0aiIILQAAIQsgBUFAayAILQABEAEgAiALOgAAIAYgBUEoaiAHEAJBAXRqIggtAAAhCyAFQShqIAgtAAEQASAEIAs6AAAgBiAFQRBqIAcQAkEBdGoiCC0AACELIAVBEGogCC0AARABIAMgCzoAACAGIAVB2ABqIAcQAkEBdGoiCC0AACELIAVB2ABqIAgtAAEQASAAIAs6AAEgBiAFQUBrIAcQAkEBdGoiCC0AACELIAVBQGsgCC0AARABIAIgCzoAASAGIAVBKGogBxACQQF0aiIILQAAIQsgBUEoaiAILQABEAEgBCALOgABIAYgBUEQaiAHEAJBAXRqIggtAAAhCyAFQRBqIAgtAAEQASADIAs6AAEgA0ECaiEDIARBAmohBCACQQJqIQIgAEECaiEAIAkgBUHYAGoQDUVxIAVBQGsQDUVxIAVBKGoQDUVxIAVBEGoQDUVxIQkMAQsLIAQgDUsgAiAMS3INAEFsIQkgACAKSw0BIApBfWohCQNAIAVB2ABqEAQgACAJT3JFBEAgBiAFQdgAaiAHEAJBAXRqIggtAAAhCyAFQdgAaiAILQABEAEgACALOgAAIAYgBUHYAGogBxACQQF0aiIILQAAIQsgBUHYAGogCC0AARABIAAgCzoAASAAQQJqIQAMAQsLA0AgBUHYAGoQBCAAIApPckUEQCAGIAVB2ABqIAcQAkEBdGoiCS0AACEIIAVB2ABqIAktAAEQASAAIAg6AAAgAEEBaiEADAELCwNAIAAgCkkEQCAGIAVB2ABqIAcQAkEBdGoiCS0AACEIIAVB2ABqIAktAAEQASAAIAg6AAAgAEEBaiEADAELCyAMQX1qIQADQCAFQUBrEAQgAiAAT3JFBEAgBiAFQUBrIAcQAkEBdGoiCi0AACEJIAVBQGsgCi0AARABIAIgCToAACAGIAVBQGsgBxACQQF0aiIKLQAAIQkgBUFAayAKLQABEAEgAiAJOgABIAJBAmohAgwBCwsDQCAFQUBrEAQgAiAMT3JFBEAgBiAFQUBrIAcQAkEBdGoiAC0AACEKIAVBQGsgAC0AARABIAIgCjoAACACQQFqIQIMAQsLA0AgAiAMSQRAIAYgBUFAayAHEAJBAXRqIgAtAAAhCiAFQUBrIAAtAAEQASACIAo6AAAgAkEBaiECDAELCyANQX1qIQADQCAFQShqEAQgBCAAT3JFBEAgBiAFQShqIAcQAkEBdGoiAi0AACEKIAVBKGogAi0AARABIAQgCjoAACAGIAVBKGogBxACQQF0aiICLQAAIQogBUEoaiACLQABEAEgBCAKOgABIARBAmohBAwBCwsDQCAFQShqEAQgBCANT3JFBEAgBiAFQShqIAcQAkEBdGoiAC0AACECIAVBKGogAC0AARABIAQgAjoAACAEQQFqIQQMAQsLA0AgBCANSQRAIAYgBUEoaiAHEAJBAXRqIgAtAAAhAiAFQShqIAAtAAEQASAEIAI6AAAgBEEBaiEEDAELCwNAIAVBEGoQBCADIA9PckUEQCAGIAVBEGogBxACQQF0aiIALQAAIQIgBUEQaiAALQABEAEgAyACOgAAIAYgBUEQaiAHEAJBAXRqIgAtAAAhAiAFQRBqIAAtAAEQASADIAI6AAEgA0ECaiEDDAELCwNAIAVBEGoQBCADIA5PckUEQCAGIAVBEGogBxACQQF0aiIALQAAIQIgBUEQaiAALQABEAEgAyACOgAAIANBAWohAwwBCwsDQCADIA5JBEAgBiAFQRBqIAcQAkEBdGoiAC0AACECIAVBEGogAC0AARABIAMgAjoAACADQQFqIQMMAQsLIAFBbCAFQdgAahAKIAVBQGsQCnEgBUEoahAKcSAFQRBqEApxGyEJDAELQWwhCQsgBUHwAGokACAJC8oCAQR/IwBBIGsiBSQAIAUgBBAOIAUtAAIhByAFQQhqIAIgAxAGIgIQA0UEQCAEQQRqIQIgACABaiIDQX1qIQQDQCAFQQhqEAQgACAET3JFBEAgAiAFQQhqIAcQAkEBdGoiBi0AACEIIAVBCGogBi0AARABIAAgCDoAACACIAVBCGogBxACQQF0aiIGLQAAIQggBUEIaiAGLQABEAEgACAIOgABIABBAmohAAwBCwsDQCAFQQhqEAQgACADT3JFBEAgAiAFQQhqIAcQAkEBdGoiBC0AACEGIAVBCGogBC0AARABIAAgBjoAACAAQQFqIQAMAQsLA0AgACADT0UEQCACIAVBCGogBxACQQF0aiIELQAAIQYgBUEIaiAELQABEAEgACAGOgAAIABBAWohAAwBCwsgAUFsIAVBCGoQChshAgsgBUEgaiQAIAILtgMBCX8jAEEQayIGJAAgBkEANgIMIAZBADYCCEFUIQQCQAJAIANBQGsiDCADIAZBCGogBkEMaiABIAIQMSICEAMNACAGQQRqIAAQDiAGKAIMIgcgBi0ABEEBaksNASAAQQRqIQogBkEAOgAFIAYgBzoABiAAIAYoAgQ2AgAgB0EBaiEJQQEhBANAIAQgCUkEQCADIARBAnRqIgEoAgAhACABIAU2AgAgACAEQX9qdCAFaiEFIARBAWohBAwBCwsgB0EBaiEHQQAhBSAGKAIIIQkDQCAFIAlGDQEgAyAFIAxqLQAAIgRBAnRqIgBBASAEdEEBdSILIAAoAgAiAWoiADYCACAHIARrIQhBACEEAkAgC0EDTQRAA0AgBCALRg0CIAogASAEakEBdGoiACAIOgABIAAgBToAACAEQQFqIQQMAAALAAsDQCABIABPDQEgCiABQQF0aiIEIAg6AAEgBCAFOgAAIAQgCDoAAyAEIAU6AAIgBCAIOgAFIAQgBToABCAEIAg6AAcgBCAFOgAGIAFBBGohAQwAAAsACyAFQQFqIQUMAAALAAsgAiEECyAGQRBqJAAgBAutAQECfwJAQYQgKAIAIABHIAAoAgBBAXYiAyABa0F4aiICQXhxQQhHcgR/IAIFIAMQJ0UNASACQQhqC0EQSQ0AIAAgACgCACICQQFxIAAgAWpBD2pBeHEiASAAa0EBdHI2AgAgASAANgIEIAEgASgCAEEBcSAAIAJBAXZqIAFrIgJBAXRyNgIAQYQgIAEgAkH/////B3FqQQRqQYQgKAIAIABGGyABNgIAIAEQJQsLygIBBX8CQAJAAkAgAEEIIABBCEsbZ0EfcyAAaUEBR2oiAUEESSAAIAF2cg0AIAFBAnRB/B5qKAIAIgJFDQADQCACQXhqIgMoAgBBAXZBeGoiBSAATwRAIAIgBUEIIAVBCEsbZ0Efc0ECdEGAH2oiASgCAEYEQCABIAIoAgQ2AgALDAMLIARBHksNASAEQQFqIQQgAigCBCICDQALC0EAIQMgAUEgTw0BA0AgAUECdEGAH2ooAgAiAkUEQCABQR5LIQIgAUEBaiEBIAJFDQEMAwsLIAIgAkF4aiIDKAIAQQF2QXhqIgFBCCABQQhLG2dBH3NBAnRBgB9qIgEoAgBGBEAgASACKAIENgIACwsgAigCACIBBEAgASACKAIENgIECyACKAIEIgEEQCABIAIoAgA2AgALIAMgAygCAEEBcjYCACADIAAQNwsgAwvhCwINfwV+IwBB8ABrIgckACAHIAAoAvDhASIINgJcIAEgAmohDSAIIAAoAoDiAWohDwJAAkAgBUUEQCABIQQMAQsgACgCxOABIRAgACgCwOABIREgACgCvOABIQ4gAEEBNgKM4QFBACEIA0AgCEEDRwRAIAcgCEECdCICaiAAIAJqQazQAWooAgA2AkQgCEEBaiEIDAELC0FsIQwgB0EYaiADIAQQBhADDQEgB0EsaiAHQRhqIAAoAgAQEyAHQTRqIAdBGGogACgCCBATIAdBPGogB0EYaiAAKAIEEBMgDUFgaiESIAEhBEEAIQwDQCAHKAIwIAcoAixBA3RqKQIAIhRCEIinQf8BcSEIIAcoAkAgBygCPEEDdGopAgAiFUIQiKdB/wFxIQsgBygCOCAHKAI0QQN0aikCACIWQiCIpyEJIBVCIIghFyAUQiCIpyECAkAgFkIQiKdB/wFxIgNBAk8EQAJAIAZFIANBGUlyRQRAIAkgB0EYaiADQSAgBygCHGsiCiAKIANLGyIKEAUgAyAKayIDdGohCSAHQRhqEAQaIANFDQEgB0EYaiADEAUgCWohCQwBCyAHQRhqIAMQBSAJaiEJIAdBGGoQBBoLIAcpAkQhGCAHIAk2AkQgByAYNwNIDAELAkAgA0UEQCACBEAgBygCRCEJDAMLIAcoAkghCQwBCwJAAkAgB0EYakEBEAUgCSACRWpqIgNBA0YEQCAHKAJEQX9qIgMgA0VqIQkMAQsgA0ECdCAHaigCRCIJIAlFaiEJIANBAUYNAQsgByAHKAJINgJMCwsgByAHKAJENgJIIAcgCTYCRAsgF6chAyALBEAgB0EYaiALEAUgA2ohAwsgCCALakEUTwRAIAdBGGoQBBoLIAgEQCAHQRhqIAgQBSACaiECCyAHQRhqEAQaIAcgB0EYaiAUQhiIp0H/AXEQCCAUp0H//wNxajYCLCAHIAdBGGogFUIYiKdB/wFxEAggFadB//8DcWo2AjwgB0EYahAEGiAHIAdBGGogFkIYiKdB/wFxEAggFqdB//8DcWo2AjQgByACNgJgIAcoAlwhCiAHIAk2AmggByADNgJkAkACQAJAIAQgAiADaiILaiASSw0AIAIgCmoiEyAPSw0AIA0gBGsgC0Egak8NAQsgByAHKQNoNwMQIAcgBykDYDcDCCAEIA0gB0EIaiAHQdwAaiAPIA4gESAQEB4hCwwBCyACIARqIQggBCAKEAcgAkERTwRAIARBEGohAgNAIAIgCkEQaiIKEAcgAkEQaiICIAhJDQALCyAIIAlrIQIgByATNgJcIAkgCCAOa0sEQCAJIAggEWtLBEBBbCELDAILIBAgAiAOayICaiIKIANqIBBNBEAgCCAKIAMQDxoMAgsgCCAKQQAgAmsQDyEIIAcgAiADaiIDNgJkIAggAmshCCAOIQILIAlBEE8EQCADIAhqIQMDQCAIIAIQByACQRBqIQIgCEEQaiIIIANJDQALDAELAkAgCUEHTQRAIAggAi0AADoAACAIIAItAAE6AAEgCCACLQACOgACIAggAi0AAzoAAyAIQQRqIAIgCUECdCIDQcAeaigCAGoiAhAXIAIgA0HgHmooAgBrIQIgBygCZCEDDAELIAggAhAMCyADQQlJDQAgAyAIaiEDIAhBCGoiCCACQQhqIgJrQQ9MBEADQCAIIAIQDCACQQhqIQIgCEEIaiIIIANJDQAMAgALAAsDQCAIIAIQByACQRBqIQIgCEEQaiIIIANJDQALCyAHQRhqEAQaIAsgDCALEAMiAhshDCAEIAQgC2ogAhshBCAFQX9qIgUNAAsgDBADDQFBbCEMIAdBGGoQBEECSQ0BQQAhCANAIAhBA0cEQCAAIAhBAnQiAmpBrNABaiACIAdqKAJENgIAIAhBAWohCAwBCwsgBygCXCEIC0G6fyEMIA8gCGsiACANIARrSw0AIAQEfyAEIAggABALIABqBUEACyABayEMCyAHQfAAaiQAIAwLkRcCFn8FfiMAQdABayIHJAAgByAAKALw4QEiCDYCvAEgASACaiESIAggACgCgOIBaiETAkACQCAFRQRAIAEhAwwBCyAAKALE4AEhESAAKALA4AEhFSAAKAK84AEhDyAAQQE2AozhAUEAIQgDQCAIQQNHBEAgByAIQQJ0IgJqIAAgAmpBrNABaigCADYCVCAIQQFqIQgMAQsLIAcgETYCZCAHIA82AmAgByABIA9rNgJoQWwhECAHQShqIAMgBBAGEAMNASAFQQQgBUEESBshFyAHQTxqIAdBKGogACgCABATIAdBxABqIAdBKGogACgCCBATIAdBzABqIAdBKGogACgCBBATQQAhBCAHQeAAaiEMIAdB5ABqIQoDQCAHQShqEARBAksgBCAXTnJFBEAgBygCQCAHKAI8QQN0aikCACIdQhCIp0H/AXEhCyAHKAJQIAcoAkxBA3RqKQIAIh5CEIinQf8BcSEJIAcoAkggBygCREEDdGopAgAiH0IgiKchCCAeQiCIISAgHUIgiKchAgJAIB9CEIinQf8BcSIDQQJPBEACQCAGRSADQRlJckUEQCAIIAdBKGogA0EgIAcoAixrIg0gDSADSxsiDRAFIAMgDWsiA3RqIQggB0EoahAEGiADRQ0BIAdBKGogAxAFIAhqIQgMAQsgB0EoaiADEAUgCGohCCAHQShqEAQaCyAHKQJUISEgByAINgJUIAcgITcDWAwBCwJAIANFBEAgAgRAIAcoAlQhCAwDCyAHKAJYIQgMAQsCQAJAIAdBKGpBARAFIAggAkVqaiIDQQNGBEAgBygCVEF/aiIDIANFaiEIDAELIANBAnQgB2ooAlQiCCAIRWohCCADQQFGDQELIAcgBygCWDYCXAsLIAcgBygCVDYCWCAHIAg2AlQLICCnIQMgCQRAIAdBKGogCRAFIANqIQMLIAkgC2pBFE8EQCAHQShqEAQaCyALBEAgB0EoaiALEAUgAmohAgsgB0EoahAEGiAHIAcoAmggAmoiCSADajYCaCAKIAwgCCAJSxsoAgAhDSAHIAdBKGogHUIYiKdB/wFxEAggHadB//8DcWo2AjwgByAHQShqIB5CGIinQf8BcRAIIB6nQf//A3FqNgJMIAdBKGoQBBogB0EoaiAfQhiIp0H/AXEQCCEOIAdB8ABqIARBBHRqIgsgCSANaiAIazYCDCALIAg2AgggCyADNgIEIAsgAjYCACAHIA4gH6dB//8DcWo2AkQgBEEBaiEEDAELCyAEIBdIDQEgEkFgaiEYIAdB4ABqIRogB0HkAGohGyABIQMDQCAHQShqEARBAksgBCAFTnJFBEAgBygCQCAHKAI8QQN0aikCACIdQhCIp0H/AXEhCyAHKAJQIAcoAkxBA3RqKQIAIh5CEIinQf8BcSEIIAcoAkggBygCREEDdGopAgAiH0IgiKchCSAeQiCIISAgHUIgiKchDAJAIB9CEIinQf8BcSICQQJPBEACQCAGRSACQRlJckUEQCAJIAdBKGogAkEgIAcoAixrIgogCiACSxsiChAFIAIgCmsiAnRqIQkgB0EoahAEGiACRQ0BIAdBKGogAhAFIAlqIQkMAQsgB0EoaiACEAUgCWohCSAHQShqEAQaCyAHKQJUISEgByAJNgJUIAcgITcDWAwBCwJAIAJFBEAgDARAIAcoAlQhCQwDCyAHKAJYIQkMAQsCQAJAIAdBKGpBARAFIAkgDEVqaiICQQNGBEAgBygCVEF/aiICIAJFaiEJDAELIAJBAnQgB2ooAlQiCSAJRWohCSACQQFGDQELIAcgBygCWDYCXAsLIAcgBygCVDYCWCAHIAk2AlQLICCnIRQgCARAIAdBKGogCBAFIBRqIRQLIAggC2pBFE8EQCAHQShqEAQaCyALBEAgB0EoaiALEAUgDGohDAsgB0EoahAEGiAHIAcoAmggDGoiGSAUajYCaCAbIBogCSAZSxsoAgAhHCAHIAdBKGogHUIYiKdB/wFxEAggHadB//8DcWo2AjwgByAHQShqIB5CGIinQf8BcRAIIB6nQf//A3FqNgJMIAdBKGoQBBogByAHQShqIB9CGIinQf8BcRAIIB+nQf//A3FqNgJEIAcgB0HwAGogBEEDcUEEdGoiDSkDCCIdNwPIASAHIA0pAwAiHjcDwAECQAJAAkAgBygCvAEiDiAepyICaiIWIBNLDQAgAyAHKALEASIKIAJqIgtqIBhLDQAgEiADayALQSBqTw0BCyAHIAcpA8gBNwMQIAcgBykDwAE3AwggAyASIAdBCGogB0G8AWogEyAPIBUgERAeIQsMAQsgAiADaiEIIAMgDhAHIAJBEU8EQCADQRBqIQIDQCACIA5BEGoiDhAHIAJBEGoiAiAISQ0ACwsgCCAdpyIOayECIAcgFjYCvAEgDiAIIA9rSwRAIA4gCCAVa0sEQEFsIQsMAgsgESACIA9rIgJqIhYgCmogEU0EQCAIIBYgChAPGgwCCyAIIBZBACACaxAPIQggByACIApqIgo2AsQBIAggAmshCCAPIQILIA5BEE8EQCAIIApqIQoDQCAIIAIQByACQRBqIQIgCEEQaiIIIApJDQALDAELAkAgDkEHTQRAIAggAi0AADoAACAIIAItAAE6AAEgCCACLQACOgACIAggAi0AAzoAAyAIQQRqIAIgDkECdCIKQcAeaigCAGoiAhAXIAIgCkHgHmooAgBrIQIgBygCxAEhCgwBCyAIIAIQDAsgCkEJSQ0AIAggCmohCiAIQQhqIgggAkEIaiICa0EPTARAA0AgCCACEAwgAkEIaiECIAhBCGoiCCAKSQ0ADAIACwALA0AgCCACEAcgAkEQaiECIAhBEGoiCCAKSQ0ACwsgCxADBEAgCyEQDAQFIA0gDDYCACANIBkgHGogCWs2AgwgDSAJNgIIIA0gFDYCBCAEQQFqIQQgAyALaiEDDAILAAsLIAQgBUgNASAEIBdrIQtBACEEA0AgCyAFSARAIAcgB0HwAGogC0EDcUEEdGoiAikDCCIdNwPIASAHIAIpAwAiHjcDwAECQAJAAkAgBygCvAEiDCAepyICaiIKIBNLDQAgAyAHKALEASIJIAJqIhBqIBhLDQAgEiADayAQQSBqTw0BCyAHIAcpA8gBNwMgIAcgBykDwAE3AxggAyASIAdBGGogB0G8AWogEyAPIBUgERAeIRAMAQsgAiADaiEIIAMgDBAHIAJBEU8EQCADQRBqIQIDQCACIAxBEGoiDBAHIAJBEGoiAiAISQ0ACwsgCCAdpyIGayECIAcgCjYCvAEgBiAIIA9rSwRAIAYgCCAVa0sEQEFsIRAMAgsgESACIA9rIgJqIgwgCWogEU0EQCAIIAwgCRAPGgwCCyAIIAxBACACaxAPIQggByACIAlqIgk2AsQBIAggAmshCCAPIQILIAZBEE8EQCAIIAlqIQYDQCAIIAIQByACQRBqIQIgCEEQaiIIIAZJDQALDAELAkAgBkEHTQRAIAggAi0AADoAACAIIAItAAE6AAEgCCACLQACOgACIAggAi0AAzoAAyAIQQRqIAIgBkECdCIGQcAeaigCAGoiAhAXIAIgBkHgHmooAgBrIQIgBygCxAEhCQwBCyAIIAIQDAsgCUEJSQ0AIAggCWohBiAIQQhqIgggAkEIaiICa0EPTARAA0AgCCACEAwgAkEIaiECIAhBCGoiCCAGSQ0ADAIACwALA0AgCCACEAcgAkEQaiECIAhBEGoiCCAGSQ0ACwsgEBADDQMgC0EBaiELIAMgEGohAwwBCwsDQCAEQQNHBEAgACAEQQJ0IgJqQazQAWogAiAHaigCVDYCACAEQQFqIQQMAQsLIAcoArwBIQgLQbp/IRAgEyAIayIAIBIgA2tLDQAgAwR/IAMgCCAAEAsgAGoFQQALIAFrIRALIAdB0AFqJAAgEAslACAAQgA3AgAgAEEAOwEIIABBADoACyAAIAE2AgwgACACOgAKC7QFAQN/IwBBMGsiBCQAIABB/wFqIgVBfWohBgJAIAMvAQIEQCAEQRhqIAEgAhAGIgIQAw0BIARBEGogBEEYaiADEBwgBEEIaiAEQRhqIAMQHCAAIQMDQAJAIARBGGoQBCADIAZPckUEQCADIARBEGogBEEYahASOgAAIAMgBEEIaiAEQRhqEBI6AAEgBEEYahAERQ0BIANBAmohAwsgBUF+aiEFAn8DQEG6fyECIAMiASAFSw0FIAEgBEEQaiAEQRhqEBI6AAAgAUEBaiEDIARBGGoQBEEDRgRAQQIhAiAEQQhqDAILIAMgBUsNBSABIARBCGogBEEYahASOgABIAFBAmohA0EDIQIgBEEYahAEQQNHDQALIARBEGoLIQUgAyAFIARBGGoQEjoAACABIAJqIABrIQIMAwsgAyAEQRBqIARBGGoQEjoAAiADIARBCGogBEEYahASOgADIANBBGohAwwAAAsACyAEQRhqIAEgAhAGIgIQAw0AIARBEGogBEEYaiADEBwgBEEIaiAEQRhqIAMQHCAAIQMDQAJAIARBGGoQBCADIAZPckUEQCADIARBEGogBEEYahAROgAAIAMgBEEIaiAEQRhqEBE6AAEgBEEYahAERQ0BIANBAmohAwsgBUF+aiEFAn8DQEG6fyECIAMiASAFSw0EIAEgBEEQaiAEQRhqEBE6AAAgAUEBaiEDIARBGGoQBEEDRgRAQQIhAiAEQQhqDAILIAMgBUsNBCABIARBCGogBEEYahAROgABIAFBAmohA0EDIQIgBEEYahAEQQNHDQALIARBEGoLIQUgAyAFIARBGGoQEToAACABIAJqIABrIQIMAgsgAyAEQRBqIARBGGoQEToAAiADIARBCGogBEEYahAROgADIANBBGohAwwAAAsACyAEQTBqJAAgAgtpAQF/An8CQAJAIAJBB00NACABKAAAQbfIwuF+Rw0AIAAgASgABDYCmOIBQWIgAEEQaiABIAIQPiIDEAMNAhogAEKBgICAEDcDiOEBIAAgASADaiACIANrECoMAQsgACABIAIQKgtBAAsLrQMBBn8jAEGAAWsiAyQAQWIhCAJAIAJBCUkNACAAQZjQAGogAUEIaiIEIAJBeGogAEGY0AAQMyIFEAMiBg0AIANBHzYCfCADIANB/ABqIANB+ABqIAQgBCAFaiAGGyIEIAEgAmoiAiAEaxAVIgUQAw0AIAMoAnwiBkEfSw0AIAMoAngiB0EJTw0AIABBiCBqIAMgBkGAC0GADCAHEBggA0E0NgJ8IAMgA0H8AGogA0H4AGogBCAFaiIEIAIgBGsQFSIFEAMNACADKAJ8IgZBNEsNACADKAJ4IgdBCk8NACAAQZAwaiADIAZBgA1B4A4gBxAYIANBIzYCfCADIANB/ABqIANB+ABqIAQgBWoiBCACIARrEBUiBRADDQAgAygCfCIGQSNLDQAgAygCeCIHQQpPDQAgACADIAZBwBBB0BEgBxAYIAQgBWoiBEEMaiIFIAJLDQAgAiAFayEFQQAhAgNAIAJBA0cEQCAEKAAAIgZBf2ogBU8NAiAAIAJBAnRqQZzQAWogBjYCACACQQFqIQIgBEEEaiEEDAELCyAEIAFrIQgLIANBgAFqJAAgCAtGAQN/IABBCGohAyAAKAIEIQJBACEAA0AgACACdkUEQCABIAMgAEEDdGotAAJBFktqIQEgAEEBaiEADAELCyABQQggAmt0C4YDAQV/Qbh/IQcCQCADRQ0AIAItAAAiBEUEQCABQQA2AgBBAUG4fyADQQFGGw8LAn8gAkEBaiIFIARBGHRBGHUiBkF/Sg0AGiAGQX9GBEAgA0EDSA0CIAUvAABBgP4BaiEEIAJBA2oMAQsgA0ECSA0BIAItAAEgBEEIdHJBgIB+aiEEIAJBAmoLIQUgASAENgIAIAVBAWoiASACIANqIgNLDQBBbCEHIABBEGogACAFLQAAIgVBBnZBI0EJIAEgAyABa0HAEEHQEUHwEiAAKAKM4QEgACgCnOIBIAQQHyIGEAMiCA0AIABBmCBqIABBCGogBUEEdkEDcUEfQQggASABIAZqIAgbIgEgAyABa0GAC0GADEGAFyAAKAKM4QEgACgCnOIBIAQQHyIGEAMiCA0AIABBoDBqIABBBGogBUECdkEDcUE0QQkgASABIAZqIAgbIgEgAyABa0GADUHgDkGQGSAAKAKM4QEgACgCnOIBIAQQHyIAEAMNACAAIAFqIAJrIQcLIAcLrQMBCn8jAEGABGsiCCQAAn9BUiACQf8BSw0AGkFUIANBDEsNABogAkEBaiELIABBBGohCUGAgAQgA0F/anRBEHUhCkEAIQJBASEEQQEgA3QiB0F/aiIMIQUDQCACIAtGRQRAAkAgASACQQF0Ig1qLwEAIgZB//8DRgRAIAkgBUECdGogAjoAAiAFQX9qIQVBASEGDAELIARBACAKIAZBEHRBEHVKGyEECyAIIA1qIAY7AQAgAkEBaiECDAELCyAAIAQ7AQIgACADOwEAIAdBA3YgB0EBdmpBA2ohBkEAIQRBACECA0AgBCALRkUEQCABIARBAXRqLgEAIQpBACEAA0AgACAKTkUEQCAJIAJBAnRqIAQ6AAIDQCACIAZqIAxxIgIgBUsNAAsgAEEBaiEADAELCyAEQQFqIQQMAQsLQX8gAg0AGkEAIQIDfyACIAdGBH9BAAUgCCAJIAJBAnRqIgAtAAJBAXRqIgEgAS8BACIBQQFqOwEAIAAgAyABEBRrIgU6AAMgACABIAVB/wFxdCAHazsBACACQQFqIQIMAQsLCyEFIAhBgARqJAAgBQvjBgEIf0FsIQcCQCACQQNJDQACQAJAAkACQCABLQAAIgNBA3EiCUEBaw4DAwEAAgsgACgCiOEBDQBBYg8LIAJBBUkNAkEDIQYgASgAACEFAn8CQAJAIANBAnZBA3EiCEF+aiIEQQFNBEAgBEEBaw0BDAILIAVBDnZB/wdxIQQgBUEEdkH/B3EhAyAIRQwCCyAFQRJ2IQRBBCEGIAVBBHZB//8AcSEDQQAMAQsgBUEEdkH//w9xIgNBgIAISw0DIAEtAARBCnQgBUEWdnIhBEEFIQZBAAshBSAEIAZqIgogAksNAgJAIANBgQZJDQAgACgCnOIBRQ0AQQAhAgNAIAJBg4ABSw0BIAJBQGshAgwAAAsACwJ/IAlBA0YEQCABIAZqIQEgAEHw4gFqIQIgACgCDCEGIAUEQCACIAMgASAEIAYQXwwCCyACIAMgASAEIAYQXQwBCyAAQbjQAWohAiABIAZqIQEgAEHw4gFqIQYgAEGo0ABqIQggBQRAIAggBiADIAEgBCACEF4MAQsgCCAGIAMgASAEIAIQXAsQAw0CIAAgAzYCgOIBIABBATYCiOEBIAAgAEHw4gFqNgLw4QEgCUECRgRAIAAgAEGo0ABqNgIMCyAAIANqIgBBiOMBakIANwAAIABBgOMBakIANwAAIABB+OIBakIANwAAIABB8OIBakIANwAAIAoPCwJ/AkACQAJAIANBAnZBA3FBf2oiBEECSw0AIARBAWsOAgACAQtBASEEIANBA3YMAgtBAiEEIAEvAABBBHYMAQtBAyEEIAEQIUEEdgsiAyAEaiIFQSBqIAJLBEAgBSACSw0CIABB8OIBaiABIARqIAMQCyEBIAAgAzYCgOIBIAAgATYC8OEBIAEgA2oiAEIANwAYIABCADcAECAAQgA3AAggAEIANwAAIAUPCyAAIAM2AoDiASAAIAEgBGo2AvDhASAFDwsCfwJAAkACQCADQQJ2QQNxQX9qIgRBAksNACAEQQFrDgIAAgELQQEhByADQQN2DAILQQIhByABLwAAQQR2DAELIAJBBEkgARAhIgJBj4CAAUtyDQFBAyEHIAJBBHYLIQIgAEHw4gFqIAEgB2otAAAgAkEgahAQIQEgACACNgKA4gEgACABNgLw4QEgB0EBaiEHCyAHC0sAIABC+erQ0OfJoeThADcDICAAQgA3AxggAELP1tO+0ser2UI3AxAgAELW64Lu6v2J9eAANwMIIABCADcDACAAQShqQQBBKBAQGgviAgICfwV+IABBKGoiASAAKAJIaiECAn4gACkDACIDQiBaBEAgACkDECIEQgeJIAApAwgiBUIBiXwgACkDGCIGQgyJfCAAKQMgIgdCEol8IAUQGSAEEBkgBhAZIAcQGQwBCyAAKQMYQsXP2bLx5brqJ3wLIAN8IQMDQCABQQhqIgAgAk0EQEIAIAEpAAAQCSADhUIbiUKHla+vmLbem55/fkLj3MqV/M7y9YV/fCEDIAAhAQwBCwsCQCABQQRqIgAgAksEQCABIQAMAQsgASgAAK1Ch5Wvr5i23puef34gA4VCF4lCz9bTvtLHq9lCfkL5893xmfaZqxZ8IQMLA0AgACACSQRAIAAxAABCxc/ZsvHluuonfiADhUILiUKHla+vmLbem55/fiEDIABBAWohAAwBCwsgA0IhiCADhULP1tO+0ser2UJ+IgNCHYggA4VC+fPd8Zn2masWfiIDQiCIIAOFC+8CAgJ/BH4gACAAKQMAIAKtfDcDAAJAAkAgACgCSCIDIAJqIgRBH00EQCABRQ0BIAAgA2pBKGogASACECAgACgCSCACaiEEDAELIAEgAmohAgJ/IAMEQCAAQShqIgQgA2ogAUEgIANrECAgACAAKQMIIAQpAAAQCTcDCCAAIAApAxAgACkAMBAJNwMQIAAgACkDGCAAKQA4EAk3AxggACAAKQMgIABBQGspAAAQCTcDICAAKAJIIQMgAEEANgJIIAEgA2tBIGohAQsgAUEgaiACTQsEQCACQWBqIQMgACkDICEFIAApAxghBiAAKQMQIQcgACkDCCEIA0AgCCABKQAAEAkhCCAHIAEpAAgQCSEHIAYgASkAEBAJIQYgBSABKQAYEAkhBSABQSBqIgEgA00NAAsgACAFNwMgIAAgBjcDGCAAIAc3AxAgACAINwMICyABIAJPDQEgAEEoaiABIAIgAWsiBBAgCyAAIAQ2AkgLCy8BAX8gAEUEQEG2f0EAIAMbDwtBun8hBCADIAFNBH8gACACIAMQEBogAwVBun8LCy8BAX8gAEUEQEG2f0EAIAMbDwtBun8hBCADIAFNBH8gACACIAMQCxogAwVBun8LC6gCAQZ/IwBBEGsiByQAIABB2OABaikDAEKAgIAQViEIQbh/IQUCQCAEQf//B0sNACAAIAMgBBBCIgUQAyIGDQAgACgCnOIBIQkgACAHQQxqIAMgAyAFaiAGGyIKIARBACAFIAYbayIGEEAiAxADBEAgAyEFDAELIAcoAgwhBCABRQRAQbp/IQUgBEEASg0BCyAGIANrIQUgAyAKaiEDAkAgCQRAIABBADYCnOIBDAELAkACQAJAIARBBUgNACAAQdjgAWopAwBCgICACFgNAAwBCyAAQQA2ApziAQwBCyAAKAIIED8hBiAAQQA2ApziASAGQRRPDQELIAAgASACIAMgBSAEIAgQOSEFDAELIAAgASACIAMgBSAEIAgQOiEFCyAHQRBqJAAgBQtnACAAQdDgAWogASACIAAoAuzhARAuIgEQAwRAIAEPC0G4fyECAkAgAQ0AIABB7OABaigCACIBBEBBYCECIAAoApjiASABRw0BC0EAIQIgAEHw4AFqKAIARQ0AIABBkOEBahBDCyACCycBAX8QVyIERQRAQUAPCyAEIAAgASACIAMgBBBLEE8hACAEEFYgAAs/AQF/AkACQAJAIAAoAqDiAUEBaiIBQQJLDQAgAUEBaw4CAAECCyAAEDBBAA8LIABBADYCoOIBCyAAKAKU4gELvAMCB38BfiMAQRBrIgkkAEG4fyEGAkAgBCgCACIIQQVBCSAAKALs4QEiBRtJDQAgAygCACIHQQFBBSAFGyAFEC8iBRADBEAgBSEGDAELIAggBUEDakkNACAAIAcgBRBJIgYQAw0AIAEgAmohCiAAQZDhAWohCyAIIAVrIQIgBSAHaiEHIAEhBQNAIAcgAiAJECwiBhADDQEgAkF9aiICIAZJBEBBuH8hBgwCCyAJKAIAIghBAksEQEFsIQYMAgsgB0EDaiEHAn8CQAJAAkAgCEEBaw4CAgABCyAAIAUgCiAFayAHIAYQSAwCCyAFIAogBWsgByAGEEcMAQsgBSAKIAVrIActAAAgCSgCCBBGCyIIEAMEQCAIIQYMAgsgACgC8OABBEAgCyAFIAgQRQsgAiAGayECIAYgB2ohByAFIAhqIQUgCSgCBEUNAAsgACkD0OABIgxCf1IEQEFsIQYgDCAFIAFrrFINAQsgACgC8OABBEBBaiEGIAJBBEkNASALEEQhDCAHKAAAIAynRw0BIAdBBGohByACQXxqIQILIAMgBzYCACAEIAI2AgAgBSABayEGCyAJQRBqJAAgBgsuACAAECsCf0EAQQAQAw0AGiABRSACRXJFBEBBYiAAIAEgAhA9EAMNARoLQQALCzcAIAEEQCAAIAAoAsTgASABKAIEIAEoAghqRzYCnOIBCyAAECtBABADIAFFckUEQCAAIAEQWwsL0QIBB38jAEEQayIGJAAgBiAENgIIIAYgAzYCDCAFBEAgBSgCBCEKIAUoAgghCQsgASEIAkACQANAIAAoAuzhARAWIQsCQANAIAQgC0kNASADKAAAQXBxQdDUtMIBRgRAIAMgBBAiIgcQAw0EIAQgB2shBCADIAdqIQMMAQsLIAYgAzYCDCAGIAQ2AggCQCAFBEAgACAFEE5BACEHQQAQA0UNAQwFCyAAIAogCRBNIgcQAw0ECyAAIAgQUCAMQQFHQQAgACAIIAIgBkEMaiAGQQhqEEwiByIDa0EAIAMQAxtBCkdyRQRAQbh/IQcMBAsgBxADDQMgAiAHayECIAcgCGohCEEBIQwgBigCDCEDIAYoAgghBAwBCwsgBiADNgIMIAYgBDYCCEG4fyEHIAQNASAIIAFrIQcMAQsgBiADNgIMIAYgBDYCCAsgBkEQaiQAIAcLRgECfyABIAAoArjgASICRwRAIAAgAjYCxOABIAAgATYCuOABIAAoArzgASEDIAAgATYCvOABIAAgASADIAJrajYCwOABCwutAgIEfwF+IwBBQGoiBCQAAkACQCACQQhJDQAgASgAAEFwcUHQ1LTCAUcNACABIAIQIiEBIABCADcDCCAAQQA2AgQgACABNgIADAELIARBGGogASACEC0iAxADBEAgACADEBoMAQsgAwRAIABBuH8QGgwBCyACIAQoAjAiA2shAiABIANqIQMDQAJAIAAgAyACIARBCGoQLCIFEAMEfyAFBSACIAVBA2oiBU8NAUG4fwsQGgwCCyAGQQFqIQYgAiAFayECIAMgBWohAyAEKAIMRQ0ACyAEKAI4BEAgAkEDTQRAIABBuH8QGgwCCyADQQRqIQMLIAQoAighAiAEKQMYIQcgAEEANgIEIAAgAyABazYCACAAIAIgBmytIAcgB0J/URs3AwgLIARBQGskAAslAQF/IwBBEGsiAiQAIAIgACABEFEgAigCACEAIAJBEGokACAAC30BBH8jAEGQBGsiBCQAIARB/wE2AggCQCAEQRBqIARBCGogBEEMaiABIAIQFSIGEAMEQCAGIQUMAQtBVCEFIAQoAgwiB0EGSw0AIAMgBEEQaiAEKAIIIAcQQSIFEAMNACAAIAEgBmogAiAGayADEDwhBQsgBEGQBGokACAFC4cBAgJ/An5BABAWIQMCQANAIAEgA08EQAJAIAAoAABBcHFB0NS0wgFGBEAgACABECIiAhADRQ0BQn4PCyAAIAEQVSIEQn1WDQMgBCAFfCIFIARUIQJCfiEEIAINAyAAIAEQUiICEAMNAwsgASACayEBIAAgAmohAAwBCwtCfiAFIAEbIQQLIAQLPwIBfwF+IwBBMGsiAiQAAn5CfiACQQhqIAAgARAtDQAaQgAgAigCHEEBRg0AGiACKQMICyEDIAJBMGokACADC40BAQJ/IwBBMGsiASQAAkAgAEUNACAAKAKI4gENACABIABB/OEBaigCADYCKCABIAApAvThATcDICAAEDAgACgCqOIBIQIgASABKAIoNgIYIAEgASkDIDcDECACIAFBEGoQGyAAQQA2AqjiASABIAEoAig2AgggASABKQMgNwMAIAAgARAbCyABQTBqJAALKgECfyMAQRBrIgAkACAAQQA2AgggAEIANwMAIAAQWCEBIABBEGokACABC4cBAQN/IwBBEGsiAiQAAkAgACgCAEUgACgCBEVzDQAgAiAAKAIINgIIIAIgACkCADcDAAJ/IAIoAgAiAQRAIAIoAghBqOMJIAERBQAMAQtBqOMJECgLIgFFDQAgASAAKQIANwL04QEgAUH84QFqIAAoAgg2AgAgARBZIAEhAwsgAkEQaiQAIAMLywEBAn8jAEEgayIBJAAgAEGBgIDAADYCtOIBIABBADYCiOIBIABBADYC7OEBIABCADcDkOIBIABBADYCpOMJIABBADYC3OIBIABCADcCzOIBIABBADYCvOIBIABBADYCxOABIABCADcCnOIBIABBpOIBakIANwIAIABBrOIBakEANgIAIAFCADcCECABQgA3AhggASABKQMYNwMIIAEgASkDEDcDACABKAIIQQh2QQFxIQIgAEEANgLg4gEgACACNgKM4gEgAUEgaiQAC3YBA38jAEEwayIBJAAgAARAIAEgAEHE0AFqIgIoAgA2AiggASAAKQK80AE3AyAgACgCACEDIAEgAigCADYCGCABIAApArzQATcDECADIAFBEGoQGyABIAEoAig2AgggASABKQMgNwMAIAAgARAbCyABQTBqJAALzAEBAX8gACABKAK00AE2ApjiASAAIAEoAgQiAjYCwOABIAAgAjYCvOABIAAgAiABKAIIaiICNgK44AEgACACNgLE4AEgASgCuNABBEAgAEKBgICAEDcDiOEBIAAgAUGk0ABqNgIMIAAgAUGUIGo2AgggACABQZwwajYCBCAAIAFBDGo2AgAgAEGs0AFqIAFBqNABaigCADYCACAAQbDQAWogAUGs0AFqKAIANgIAIABBtNABaiABQbDQAWooAgA2AgAPCyAAQgA3A4jhAQs7ACACRQRAQbp/DwsgBEUEQEFsDwsgAiAEEGAEQCAAIAEgAiADIAQgBRBhDwsgACABIAIgAyAEIAUQZQtGAQF/IwBBEGsiBSQAIAVBCGogBBAOAn8gBS0ACQRAIAAgASACIAMgBBAyDAELIAAgASACIAMgBBA0CyEAIAVBEGokACAACzQAIAAgAyAEIAUQNiIFEAMEQCAFDwsgBSAESQR/IAEgAiADIAVqIAQgBWsgABA1BUG4fwsLRgEBfyMAQRBrIgUkACAFQQhqIAQQDgJ/IAUtAAkEQCAAIAEgAiADIAQQYgwBCyAAIAEgAiADIAQQNQshACAFQRBqJAAgAAtZAQF/QQ8hAiABIABJBEAgAUEEdCAAbiECCyAAQQh2IgEgAkEYbCIAQYwIaigCAGwgAEGICGooAgBqIgJBA3YgAmogAEGACGooAgAgAEGECGooAgAgAWxqSQs3ACAAIAMgBCAFQYAQEDMiBRADBEAgBQ8LIAUgBEkEfyABIAIgAyAFaiAEIAVrIAAQMgVBuH8LC78DAQN/IwBBIGsiBSQAIAVBCGogAiADEAYiAhADRQRAIAAgAWoiB0F9aiEGIAUgBBAOIARBBGohAiAFLQACIQMDQEEAIAAgBkkgBUEIahAEGwRAIAAgAiAFQQhqIAMQAkECdGoiBC8BADsAACAFQQhqIAQtAAIQASAAIAQtAANqIgQgAiAFQQhqIAMQAkECdGoiAC8BADsAACAFQQhqIAAtAAIQASAEIAAtAANqIQAMAQUgB0F+aiEEA0AgBUEIahAEIAAgBEtyRQRAIAAgAiAFQQhqIAMQAkECdGoiBi8BADsAACAFQQhqIAYtAAIQASAAIAYtAANqIQAMAQsLA0AgACAES0UEQCAAIAIgBUEIaiADEAJBAnRqIgYvAQA7AAAgBUEIaiAGLQACEAEgACAGLQADaiEADAELCwJAIAAgB08NACAAIAIgBUEIaiADEAIiA0ECdGoiAC0AADoAACAALQADQQFGBEAgBUEIaiAALQACEAEMAQsgBSgCDEEfSw0AIAVBCGogAiADQQJ0ai0AAhABIAUoAgxBIUkNACAFQSA2AgwLIAFBbCAFQQhqEAobIQILCwsgBUEgaiQAIAILkgIBBH8jAEFAaiIJJAAgCSADQTQQCyEDAkAgBEECSA0AIAMgBEECdGooAgAhCSADQTxqIAgQIyADQQE6AD8gAyACOgA+QQAhBCADKAI8IQoDQCAEIAlGDQEgACAEQQJ0aiAKNgEAIARBAWohBAwAAAsAC0EAIQkDQCAGIAlGRQRAIAMgBSAJQQF0aiIKLQABIgtBAnRqIgwoAgAhBCADQTxqIAotAABBCHQgCGpB//8DcRAjIANBAjoAPyADIAcgC2siCiACajoAPiAEQQEgASAKa3RqIQogAygCPCELA0AgACAEQQJ0aiALNgEAIARBAWoiBCAKSQ0ACyAMIAo2AgAgCUEBaiEJDAELCyADQUBrJAALowIBCX8jAEHQAGsiCSQAIAlBEGogBUE0EAsaIAcgBmshDyAHIAFrIRADQAJAIAMgCkcEQEEBIAEgByACIApBAXRqIgYtAAEiDGsiCGsiC3QhDSAGLQAAIQ4gCUEQaiAMQQJ0aiIMKAIAIQYgCyAPTwRAIAAgBkECdGogCyAIIAUgCEE0bGogCCAQaiIIQQEgCEEBShsiCCACIAQgCEECdGooAgAiCEEBdGogAyAIayAHIA4QYyAGIA1qIQgMAgsgCUEMaiAOECMgCUEBOgAPIAkgCDoADiAGIA1qIQggCSgCDCELA0AgBiAITw0CIAAgBkECdGogCzYBACAGQQFqIQYMAAALAAsgCUHQAGokAA8LIAwgCDYCACAKQQFqIQoMAAALAAs0ACAAIAMgBCAFEDYiBRADBEAgBQ8LIAUgBEkEfyABIAIgAyAFaiAEIAVrIAAQNAVBuH8LCyMAIAA/AEEQdGtB//8DakEQdkAAQX9GBEBBAA8LQQAQAEEBCzsBAX8gAgRAA0AgACABIAJBgCAgAkGAIEkbIgMQCyEAIAFBgCBqIQEgAEGAIGohACACIANrIgINAAsLCwYAIAAQAwsLqBUJAEGICAsNAQAAAAEAAAACAAAAAgBBoAgLswYBAAAAAQAAAAIAAAACAAAAJgAAAIIAAAAhBQAASgAAAGcIAAAmAAAAwAEAAIAAAABJBQAASgAAAL4IAAApAAAALAIAAIAAAABJBQAASgAAAL4IAAAvAAAAygIAAIAAAACKBQAASgAAAIQJAAA1AAAAcwMAAIAAAACdBQAASgAAAKAJAAA9AAAAgQMAAIAAAADrBQAASwAAAD4KAABEAAAAngMAAIAAAABNBgAASwAAAKoKAABLAAAAswMAAIAAAADBBgAATQAAAB8NAABNAAAAUwQAAIAAAAAjCAAAUQAAAKYPAABUAAAAmQQAAIAAAABLCQAAVwAAALESAABYAAAA2gQAAIAAAABvCQAAXQAAACMUAABUAAAARQUAAIAAAABUCgAAagAAAIwUAABqAAAArwUAAIAAAAB2CQAAfAAAAE4QAAB8AAAA0gIAAIAAAABjBwAAkQAAAJAHAACSAAAAAAAAAAEAAAABAAAABQAAAA0AAAAdAAAAPQAAAH0AAAD9AAAA/QEAAP0DAAD9BwAA/Q8AAP0fAAD9PwAA/X8AAP3/AAD9/wEA/f8DAP3/BwD9/w8A/f8fAP3/PwD9/38A/f//AP3//wH9//8D/f//B/3//w/9//8f/f//P/3//38AAAAAAQAAAAIAAAADAAAABAAAAAUAAAAGAAAABwAAAAgAAAAJAAAACgAAAAsAAAAMAAAADQAAAA4AAAAPAAAAEAAAABEAAAASAAAAEwAAABQAAAAVAAAAFgAAABcAAAAYAAAAGQAAABoAAAAbAAAAHAAAAB0AAAAeAAAAHwAAAAMAAAAEAAAABQAAAAYAAAAHAAAACAAAAAkAAAAKAAAACwAAAAwAAAANAAAADgAAAA8AAAAQAAAAEQAAABIAAAATAAAAFAAAABUAAAAWAAAAFwAAABgAAAAZAAAAGgAAABsAAAAcAAAAHQAAAB4AAAAfAAAAIAAAACEAAAAiAAAAIwAAACUAAAAnAAAAKQAAACsAAAAvAAAAMwAAADsAAABDAAAAUwAAAGMAAACDAAAAAwEAAAMCAAADBAAAAwgAAAMQAAADIAAAA0AAAAOAAAADAAEAQeAPC1EBAAAAAQAAAAEAAAABAAAAAgAAAAIAAAADAAAAAwAAAAQAAAAEAAAABQAAAAcAAAAIAAAACQAAAAoAAAALAAAADAAAAA0AAAAOAAAADwAAABAAQcQQC4sBAQAAAAIAAAADAAAABAAAAAUAAAAGAAAABwAAAAgAAAAJAAAACgAAAAsAAAAMAAAADQAAAA4AAAAPAAAAEAAAABIAAAAUAAAAFgAAABgAAAAcAAAAIAAAACgAAAAwAAAAQAAAAIAAAAAAAQAAAAIAAAAEAAAACAAAABAAAAAgAAAAQAAAAIAAAAAAAQBBkBIL5gQBAAAAAQAAAAEAAAABAAAAAgAAAAIAAAADAAAAAwAAAAQAAAAGAAAABwAAAAgAAAAJAAAACgAAAAsAAAAMAAAADQAAAA4AAAAPAAAAEAAAAAEAAAAEAAAACAAAAAAAAAABAAEBBgAAAAAAAAQAAAAAEAAABAAAAAAgAAAFAQAAAAAAAAUDAAAAAAAABQQAAAAAAAAFBgAAAAAAAAUHAAAAAAAABQkAAAAAAAAFCgAAAAAAAAUMAAAAAAAABg4AAAAAAAEFEAAAAAAAAQUUAAAAAAABBRYAAAAAAAIFHAAAAAAAAwUgAAAAAAAEBTAAAAAgAAYFQAAAAAAABwWAAAAAAAAIBgABAAAAAAoGAAQAAAAADAYAEAAAIAAABAAAAAAAAAAEAQAAAAAAAAUCAAAAIAAABQQAAAAAAAAFBQAAACAAAAUHAAAAAAAABQgAAAAgAAAFCgAAAAAAAAULAAAAAAAABg0AAAAgAAEFEAAAAAAAAQUSAAAAIAABBRYAAAAAAAIFGAAAACAAAwUgAAAAAAADBSgAAAAAAAYEQAAAABAABgRAAAAAIAAHBYAAAAAAAAkGAAIAAAAACwYACAAAMAAABAAAAAAQAAAEAQAAACAAAAUCAAAAIAAABQMAAAAgAAAFBQAAACAAAAUGAAAAIAAABQgAAAAgAAAFCQAAACAAAAULAAAAIAAABQwAAAAAAAAGDwAAACAAAQUSAAAAIAABBRQAAAAgAAIFGAAAACAAAgUcAAAAIAADBSgAAAAgAAQFMAAAAAAAEAYAAAEAAAAPBgCAAAAAAA4GAEAAAAAADQYAIABBgBcLhwIBAAEBBQAAAAAAAAUAAAAAAAAGBD0AAAAAAAkF/QEAAAAADwX9fwAAAAAVBf3/HwAAAAMFBQAAAAAABwR9AAAAAAAMBf0PAAAAABIF/f8DAAAAFwX9/38AAAAFBR0AAAAAAAgE/QAAAAAADgX9PwAAAAAUBf3/DwAAAAIFAQAAABAABwR9AAAAAAALBf0HAAAAABEF/f8BAAAAFgX9/z8AAAAEBQ0AAAAQAAgE/QAAAAAADQX9HwAAAAATBf3/BwAAAAEFAQAAABAABgQ9AAAAAAAKBf0DAAAAABAF/f8AAAAAHAX9//8PAAAbBf3//wcAABoF/f//AwAAGQX9//8BAAAYBf3//wBBkBkLhgQBAAEBBgAAAAAAAAYDAAAAAAAABAQAAAAgAAAFBQAAAAAAAAUGAAAAAAAABQgAAAAAAAAFCQAAAAAAAAULAAAAAAAABg0AAAAAAAAGEAAAAAAAAAYTAAAAAAAABhYAAAAAAAAGGQAAAAAAAAYcAAAAAAAABh8AAAAAAAAGIgAAAAAAAQYlAAAAAAABBikAAAAAAAIGLwAAAAAAAwY7AAAAAAAEBlMAAAAAAAcGgwAAAAAACQYDAgAAEAAABAQAAAAAAAAEBQAAACAAAAUGAAAAAAAABQcAAAAgAAAFCQAAAAAAAAUKAAAAAAAABgwAAAAAAAAGDwAAAAAAAAYSAAAAAAAABhUAAAAAAAAGGAAAAAAAAAYbAAAAAAAABh4AAAAAAAAGIQAAAAAAAQYjAAAAAAABBicAAAAAAAIGKwAAAAAAAwYzAAAAAAAEBkMAAAAAAAUGYwAAAAAACAYDAQAAIAAABAQAAAAwAAAEBAAAABAAAAQFAAAAIAAABQcAAAAgAAAFCAAAACAAAAUKAAAAIAAABQsAAAAAAAAGDgAAAAAAAAYRAAAAAAAABhQAAAAAAAAGFwAAAAAAAAYaAAAAAAAABh0AAAAAAAAGIAAAAAAAEAYDAAEAAAAPBgOAAAAAAA4GA0AAAAAADQYDIAAAAAAMBgMQAAAAAAsGAwgAAAAACgYDBABBpB0L2QEBAAAAAwAAAAcAAAAPAAAAHwAAAD8AAAB/AAAA/wAAAP8BAAD/AwAA/wcAAP8PAAD/HwAA/z8AAP9/AAD//wAA//8BAP//AwD//wcA//8PAP//HwD//z8A//9/AP///wD///8B////A////wf///8P////H////z////9/AAAAAAEAAAACAAAABAAAAAAAAAACAAAABAAAAAgAAAAAAAAAAQAAAAIAAAABAAAABAAAAAQAAAAEAAAABAAAAAgAAAAIAAAACAAAAAcAAAAIAAAACQAAAAoAAAALAEGgIAsDwBBQ",
	e = [171, 75, 84, 88, 32, 50, 48, 187, 13, 10, 26, 10];
var n, i, s, a, r, o, l, f;
! function(e) {
	e[e.NONE = 0] = "NONE", e[e.BASISLZ = 1] = "BASISLZ", e[e.ZSTD = 2] = "ZSTD", e[e.ZLIB = 3] = "ZLIB"
}(n || (n = {})),
function(e) {
	e[e.BASICFORMAT = 0] = "BASICFORMAT"
}(i || (i = {})),
function(e) {
	e[e.UNSPECIFIED = 0] = "UNSPECIFIED", e[e.ETC1S = 163] = "ETC1S", e[e.UASTC = 166] = "UASTC"
}(s || (s = {})),
function(e) {
	e[e.UNSPECIFIED = 0] = "UNSPECIFIED", e[e.SRGB = 1] = "SRGB"
}(a || (a = {})),
function(e) {
	e[e.UNSPECIFIED = 0] = "UNSPECIFIED", e[e.LINEAR = 1] = "LINEAR", e[e.SRGB = 2] = "SRGB", e[e.ITU = 3] = "ITU", e[e.NTSC = 4] = "NTSC", e[e.SLOG = 5] = "SLOG", e[e.SLOG2 = 6] = "SLOG2"
}(r || (r = {})),
function(e) {
	e[e.ALPHA_STRAIGHT = 0] = "ALPHA_STRAIGHT", e[e.ALPHA_PREMULTIPLIED = 1] = "ALPHA_PREMULTIPLIED"
}(o || (o = {})),
function(e) {
	e[e.RGB = 0] = "RGB", e[e.RRR = 3] = "RRR", e[e.GGG = 4] = "GGG", e[e.AAA = 15] = "AAA"
}(l || (l = {})),
function(e) {
	e[e.RGB = 0] = "RGB", e[e.RGBA = 3] = "RGBA", e[e.RRR = 4] = "RRR", e[e.RRRG = 5] = "RRRG"
}(f || (f = {}));
class U {
	constructor() {
		this.vkFormat = 0, this.typeSize = 1, this.pixelWidth = 0, this.pixelHeight = 0, this.pixelDepth = 0, this.layerCount = 0, this.faceCount = 1, this.supercompressionScheme = n.NONE, this.levels = [], this.dataFormatDescriptor = [{
			vendorId: 0,
			descriptorType: i.BASICFORMAT,
			versionNumber: 2,
			descriptorBlockSize: 40,
			colorModel: s.UNSPECIFIED,
			colorPrimaries: a.SRGB,
			transferFunction: a.SRGB,
			flags: o.ALPHA_STRAIGHT,
			texelBlockDimension: {
				x: 4,
				y: 4,
				z: 1,
				w: 1
			},
			bytesPlane: [],
			samples: []
		}], this.keyValue = {}, this.globalData = null
	}
}
class c {
	constructor(e, t, i, n) {
		this._dataView = new DataView(e.buffer, e.byteOffset + t, i), this._littleEndian = n, this._offset = 0
	}
	_nextUint8() {
		const e = this._dataView.getUint8(this._offset);
		return this._offset += 1, e
	}
	_nextUint16() {
		const e = this._dataView.getUint16(this._offset, this._littleEndian);
		return this._offset += 2, e
	}
	_nextUint32() {
		const e = this._dataView.getUint32(this._offset, this._littleEndian);
		return this._offset += 4, e
	}
	_nextUint64() {
		const e = this._dataView.getUint32(this._offset, this._littleEndian) + 2 ** 32 * this._dataView.getUint32(this._offset + 4, this._littleEndian);
		return this._offset += 8, e
	}
	_skip(e) {
		return this._offset += e, this
	}
	_scan(e, t = 0) {
		const i = this._offset;
		let n = 0;
		for (; this._dataView.getUint8(this._offset) !== t && n < e;) n++, this._offset++;
		return n < e && this._offset++, new Uint8Array(this._dataView.buffer, this._dataView.byteOffset + i, n)
	}
}

function _(e) {
	return "undefined" != typeof TextDecoder ? (new TextDecoder).decode(e) : Buffer.from(e).toString("utf8")
}

function p(t) {
	const i = new Uint8Array(t.buffer, t.byteOffset, e.length);
	if (i[0] !== e[0] || i[1] !== e[1] || i[2] !== e[2] || i[3] !== e[3] || i[4] !== e[4] || i[5] !== e[5] || i[6] !== e[6] || i[7] !== e[7] || i[8] !== e[8] || i[9] !== e[9] || i[10] !== e[10] || i[11] !== e[11]) throw new Error("Missing KTX 2.0 identifier.");
	const n = new U,
		r = 17 * Uint32Array.BYTES_PER_ELEMENT,
		a = new c(t, e.length, r, !0);
	n.vkFormat = a._nextUint32(), n.typeSize = a._nextUint32(), n.pixelWidth = a._nextUint32(), n.pixelHeight = a._nextUint32(), n.pixelDepth = a._nextUint32(), n.layerCount = a._nextUint32(), n.faceCount = a._nextUint32();
	const s = a._nextUint32();
	n.supercompressionScheme = a._nextUint32();
	const o = a._nextUint32(),
		l = a._nextUint32(),
		h = a._nextUint32(),
		u = a._nextUint32(),
		d = a._nextUint64(),
		p = a._nextUint64(),
		m = new c(t, e.length + r, 3 * s * 8, !0);
	for (let e = 0; e < s; e++) n.levels.push({
		levelData: new Uint8Array(t.buffer, t.byteOffset + m._nextUint64(), m._nextUint64()),
		uncompressedByteLength: m._nextUint64()
	});
	const A = new c(t, o, l, !0),
		g = {
			vendorId: A._skip(4)._nextUint16(),
			descriptorType: A._nextUint16(),
			versionNumber: A._nextUint16(),
			descriptorBlockSize: A._nextUint16(),
			colorModel: A._nextUint8(),
			colorPrimaries: A._nextUint8(),
			transferFunction: A._nextUint8(),
			flags: A._nextUint8(),
			texelBlockDimension: {
				x: A._nextUint8() + 1,
				y: A._nextUint8() + 1,
				z: A._nextUint8() + 1,
				w: A._nextUint8() + 1
			},
			bytesPlane: [A._nextUint8(), A._nextUint8(), A._nextUint8(), A._nextUint8(), A._nextUint8(), A._nextUint8(), A._nextUint8(), A._nextUint8()],
			samples: []
		},
		f = (g.descriptorBlockSize / 4 - 6) / 4;
	for (let e = 0; e < f; e++) g.samples[e] = {
		bitOffset: A._nextUint16(),
		bitLength: A._nextUint8(),
		channelID: A._nextUint8(),
		samplePosition: [A._nextUint8(), A._nextUint8(), A._nextUint8(), A._nextUint8()],
		sampleLower: A._nextUint32(),
		sampleUpper: A._nextUint32()
	};
	n.dataFormatDescriptor.length = 0, n.dataFormatDescriptor.push(g);
	const v = new c(t, h, u, !0);
	for (; v._offset < u;) {
		const e = v._nextUint32(),
			t = v._scan(e),
			i = _(t),
			r = v._scan(e - t.byteLength);
		n.keyValue[i] = i.match(/^ktx/i) ? _(r) : r, v._offset % 4 && v._skip(4 - v._offset % 4)
	}
	if (p <= 0) return n;
	const y = new c(t, d, p, !0),
		E = y._nextUint16(),
		b = y._nextUint16(),
		x = y._nextUint32(),
		w = y._nextUint32(),
		C = y._nextUint32(),
		S = y._nextUint32(),
		I = [];
	for (let e = 0; e < s; e++) I.push({
		imageFlags: y._nextUint32(),
		rgbSliceByteOffset: y._nextUint32(),
		rgbSliceByteLength: y._nextUint32(),
		alphaSliceByteOffset: y._nextUint32(),
		alphaSliceByteLength: y._nextUint32()
	});
	const M = d + y._offset,
		T = M + x,
		B = T + w,
		L = B + C,
		R = new Uint8Array(t.buffer, t.byteOffset + M, x),
		D = new Uint8Array(t.buffer, t.byteOffset + T, w),
		P = new Uint8Array(t.buffer, t.byteOffset + B, C),
		Q = new Uint8Array(t.buffer, t.byteOffset + L, S);
	return n.globalData = {
		endpointCount: E,
		selectorCount: b,
		imageDescs: I,
		endpointsData: R,
		selectorsData: D,
		tablesData: P,
		extendedData: Q
	}, n
}
var DFDModel = {
		ETC1S: 163,
		UASTC: 166
	},
	DFDChannel = {
		ETC1S: {
			RGB: 0,
			RRR: 3,
			GGG: 4,
			AAA: 15
		},
		UASTC: {
			RGB: 0,
			RGBA: 3,
			RRR: 4,
			RRRG: 5
		}
	},
	SupercompressionScheme = {
		ZSTD: 2
	},
	Transfer = {
		SRGB: 2
	};
class KTX2Loader extends CompressedTextureLoader {
	constructor(e) {
		super(e), this.basisLoader = new BasisTextureLoader(e), this.zstd = new ZSTDDecoder, this.zstd.init(), "undefined" != typeof MSC_TRANSCODER && console.warn('THREE.KTX2Loader: Please update to latest "basis_transcoder". "msc_basis_transcoder" is no longer supported in three.js r125+.')
	}
	setTranscoderPath(e) {
		return this.basisLoader.setTranscoderPath(e), this
	}
	setWorkerLimit(e) {
		return this.basisLoader.setWorkerLimit(e), this
	}
	detectSupport(e) {
		return this.basisLoader.detectSupport(e), this
	}
	dispose() {
		return this.basisLoader.dispose(), this
	}
	load(e, t, i, n) {
		var r = this,
			a = new CompressedTexture;
		return new Promise((function(t, n) {
			new FileLoader(r.manager).setPath(r.path).setResponseType("arraybuffer").load(e, t, i, n)
		})).then((function(e) {
			r.parse(e, (function(e) {
				a.copy(e), a.needsUpdate = !0, t && t(a)
			}), n)
		})).catch(n), a
	}
	parse(e, t, i) {
		var n = this,
			r = p(new Uint8Array(e));
		if (r.pixelDepth > 0) throw new Error("THREE.KTX2Loader: Only 2D textures are currently supported.");
		if (r.layerCount > 1) throw new Error("THREE.KTX2Loader: Array textures are not currently supported.");
		if (r.faceCount > 1) throw new Error("THREE.KTX2Loader: Cube textures are not currently supported.");
		var a = KTX2Utils.getBasicDFD(r);
		return KTX2Utils.createLevels(r, this.zstd).then((function(e) {
			var t = a.colorModel === DFDModel.UASTC ? BasisTextureLoader.BasisFormat.UASTC_4x4 : BasisTextureLoader.BasisFormat.ETC1S,
				i = {
					levels: e,
					width: r.pixelWidth,
					height: r.pixelHeight,
					basisFormat: t,
					hasAlpha: KTX2Utils.getAlpha(r)
				};
			return t === BasisTextureLoader.BasisFormat.ETC1S && (i.globalData = r.globalData), n.basisLoader.parseInternalAsync(i)
		})).then((function(e) {
			e.encoding = a.transferFunction === Transfer.SRGB ? 3001 : 3e3, e.premultiplyAlpha = KTX2Utils.getPremultiplyAlpha(r), t(e)
		})).catch(i), this
	}
}
var _a$b, _b$9, KTX2Utils = {
	createLevels: async function(e, t) {
		e.supercompressionScheme === SupercompressionScheme.ZSTD && await t.init();
		for (var i = [], n = e.pixelWidth, r = e.pixelHeight, a = 0; a < e.levels.length; a++) {
			var s = Math.max(1, Math.floor(n / Math.pow(2, a))),
				o = Math.max(1, Math.floor(r / Math.pow(2, a))),
				l = e.levels[a].levelData;
			e.supercompressionScheme === SupercompressionScheme.ZSTD && (l = t.decode(l, e.levels[a].uncompressedByteLength)), i.push({
				index: a,
				width: s,
				height: o,
				data: l
			})
		}
		return i
	},
	getBasicDFD: function(e) {
		return e.dataFormatDescriptor[0]
	},
	getAlpha: function(e) {
		var t = this.getBasicDFD(e);
		return t.colorModel === DFDModel.UASTC ? (15 & t.samples[0].channelID) === DFDChannel.UASTC.RGBA : 2 === t.samples.length && (15 & t.samples[1].channelID) === DFDChannel.ETC1S.AAA
	},
	getPremultiplyAlpha: function(e) {
		return !!(1 & this.getBasicDFD(e).flags)
	}
};
const $retainerCount = Symbol("retainerCount"),
	$recentlyUsed = Symbol("recentlyUsed"),
	$evict = Symbol("evict"),
	$evictionThreshold = Symbol("evictionThreshold"),
	$cache = Symbol("cache");
class CacheEvictionPolicy {
	constructor(e, t = 5) {
		this[_a$b] = new Map, this[_b$9] = [], this[$cache] = e, this[$evictionThreshold] = t
	}
	set evictionThreshold(e) {
		this[$evictionThreshold] = e, this[$evict]()
	}
	get evictionThreshold() {
		return this[$evictionThreshold]
	}
	get cache() {
		return this[$cache]
	}
	retainerCount(e) {
		return this[$retainerCount].get(e) || 0
	}
	reset() {
		this[$retainerCount].clear(), this[$recentlyUsed] = []
	}
	retain(e) {
		this[$retainerCount].has(e) || this[$retainerCount].set(e, 0), this[$retainerCount].set(e, this[$retainerCount].get(e) + 1);
		const t = this[$recentlyUsed].indexOf(e); - 1 !== t && this[$recentlyUsed].splice(t, 1), this[$recentlyUsed].unshift(e), this[$evict]()
	}
	release(e) {
		this[$retainerCount].has(e) && this[$retainerCount].set(e, Math.max(this[$retainerCount].get(e) - 1, 0)), this[$evict]()
	} [(_a$b = $retainerCount, _b$9 = $recentlyUsed, $evict)]() {
		if (!(this[$recentlyUsed].length < this[$evictionThreshold]))
			for (let e = this[$recentlyUsed].length - 1; e >= this[$evictionThreshold]; --e) {
				const t = this[$recentlyUsed][e];
				0 === this[$retainerCount].get(t) && (this[$cache].delete(t), this[$recentlyUsed].splice(e, 1))
			}
	}
}
var _a$a, _b$8;
const loadWithLoader = (e, t, i = (() => {})) => {
		const n = e => {
			const t = e.loaded / e.total;
			i(Math.max(0, Math.min(1, isFinite(t) ? t : 1)))
		};
		return new Promise((i, r) => {
			t.load(e, i, n, r)
		})
	},
	cache = new Map,
	preloaded = new Map;
let dracoDecoderLocation;
const dracoLoader = new DRACOLoader;
let ktx2TranscoderLocation;
const ktx2Loader = new KTX2Loader,
	$loader = Symbol("loader"),
	$evictionPolicy = Symbol("evictionPolicy"),
	$GLTFInstance = Symbol("GLTFInstance");
class CachingGLTFLoader extends EventDispatcher {
	constructor(e) {
		super(), this[_b$8] = new GLTFLoader, this[$GLTFInstance] = e, this[$loader].setDRACOLoader(dracoLoader), this[$loader].setKTX2Loader(ktx2Loader)
	}
	static setDRACODecoderLocation(e) {
		dracoDecoderLocation = e, dracoLoader.setDecoderPath(e)
	}
	static getDRACODecoderLocation() {
		return dracoDecoderLocation
	}
	static setKTX2TranscoderLocation(e) {
		ktx2TranscoderLocation = e, ktx2Loader.setTranscoderPath(e)
	}
	static getKTX2TranscoderLocation() {
		return ktx2TranscoderLocation
	}
	static initializeKTX2Loader(e) {
		ktx2Loader.detectSupport(e)
	}
	static get cache() {
		return cache
	}
	static clearCache() {
		cache.forEach((e, t) => {
			this.delete(t)
		}), this[$evictionPolicy].reset()
	}
	static has(e) {
		return cache.has(e)
	}
	static async delete(e) {
		if (!this.has(e)) return;
		const t = cache.get(e);
		preloaded.delete(e), cache.delete(e);
		(await t).dispose()
	}
	static hasFinishedLoading(e) {
		return !!preloaded.get(e)
	}
	get[(_a$a = $evictionPolicy, _b$8 = $loader, $evictionPolicy)]() {
		return this.constructor[$evictionPolicy]
	}
	async preload(e, t, i = (() => {})) {
		if (this.dispatchEvent({
				type: "preload",
				element: t,
				src: e
			}), !cache.has(e)) {
			const t = loadWithLoader(e, this[$loader], e => {
					i(.8 * e)
				}),
				n = this[$GLTFInstance],
				r = t.then(e => n.prepare(e)).then(e => (i(.9), new n(e)));
			cache.set(e, r)
		}
		await cache.get(e), preloaded.set(e, !0), i && i(1)
	}
	async load(e, t, i = (() => {})) {
		await this.preload(e, t, i);
		const n = await cache.get(e),
			r = await n.clone();
		return this[$evictionPolicy].retain(e), r.dispose = (() => {
			const t = r.dispose;
			let i = !1;
			return () => {
				i || (i = !0, t.apply(r), this[$evictionPolicy].release(e))
			}
		})(), r
	}
}
var _a$9;
CachingGLTFLoader[_a$a] = new CacheEvictionPolicy(CachingGLTFLoader);
const SETTLING_TIME = 1e4,
	DECAY_MILLISECONDS = 50,
	NATURAL_FREQUENCY = .02,
	NIL_SPEED = 2e-4 * .02,
	$velocity = Symbol("velocity");
class Damper {
	constructor() {
		this[_a$9] = 0
	}
	update(e, t, i, n) {
		if (null == e || 0 === n) return t;
		if (e === t && 0 === this[$velocity]) return t;
		if (i < 0) return e;
		const r = e - t,
			a = this[$velocity] + .02 * r,
			s = r + i * a,
			o = Math.exp(-.02 * i),
			l = (a - .02 * s) * o,
			c = -.02 * (l + a * o);
		return Math.abs(l) < NIL_SPEED * Math.abs(n) && c * r >= 0 ? (this[$velocity] = 0, t) : (this[$velocity] = l, t + s * o)
	}
}
_a$9 = $velocity;
var CSS2DObject = function(e) {
	Object3D.call(this), this.element = e || document.createElement("div"), this.element.style.position = "absolute", this.addEventListener("removed", (function() {
		this.traverse((function(e) {
			e.element instanceof Element && null !== e.element.parentNode && e.element.parentNode.removeChild(e.element)
		}))
	}))
};
CSS2DObject.prototype = Object.assign(Object.create(Object3D.prototype), {
	constructor: CSS2DObject,
	copy: function(e, t) {
		return Object3D.prototype.copy.call(this, e, t), this.element = e.element.cloneNode(!0), this
	}
});
var CSS2DRenderer = function() {
	var e, t, i, n, r = this,
		a = new Vector3,
		s = new Matrix4,
		o = new Matrix4,
		l = {
			objects: new WeakMap
		},
		c = document.createElement("div");
	c.style.overflow = "hidden", this.domElement = c, this.getSize = function() {
		return {
			width: e,
			height: t
		}
	}, this.setSize = function(r, a) {
		i = (e = r) / 2, n = (t = a) / 2, c.style.width = r + "px", c.style.height = a + "px"
	};
	var h = function(e, t, s) {
			if (e instanceof CSS2DObject) {
				e.onBeforeRender(r, t, s), a.setFromMatrixPosition(e.matrixWorld), a.applyMatrix4(o);
				var d = e.element;
				d.style.transform = "translate(-50%,-50%) translate(" + (a.x * i + i) + "px," + (-a.y * n + n) + "px)", d.style.display = e.visible && a.z >= -1 && a.z <= 1 ? "" : "none";
				var p = {
					distanceToCameraSquared: u(s, e)
				};
				l.objects.set(e, p), d.parentNode !== c && c.appendChild(d), e.onAfterRender(r, t, s)
			}
			for (var m = 0, A = e.children.length; m < A; m++) h(e.children[m], t, s)
		},
		u = function() {
			var e = new Vector3,
				t = new Vector3;
			return function(i, n) {
				return e.setFromMatrixPosition(i.matrixWorld), t.setFromMatrixPosition(n.matrixWorld), e.distanceToSquared(t)
			}
		}(),
		d = function(e) {
			for (var t = function(e) {
					var t = [];
					return e.traverse((function(e) {
						e instanceof CSS2DObject && t.push(e)
					})), t
				}(e).sort((function(e, t) {
					return l.objects.get(e).distanceToCameraSquared - l.objects.get(t).distanceToCameraSquared
				})), i = t.length, n = 0, r = t.length; n < r; n++) t[n].element.style.zIndex = i - n
		};
	this.render = function(e, t) {
		!0 === e.autoUpdate && e.updateMatrixWorld(), null === t.parent && t.updateMatrixWorld(), s.copy(t.matrixWorldInverse), o.multiplyMatrices(t.projectionMatrix, s), h(e, e, t), d(e)
	}
};
const numberNode = (e, t) => ({
		type: "number",
		number: e,
		unit: t
	}),
	parseExpressions = (() => {
		const e = {};
		return t => {
			const i = t;
			if (i in e) return e[i];
			const n = [];
			let r = 0;
			for (; t;) {
				if (++r > 1e3) {
					t = "";
					break
				}
				const e = parseExpression(t),
					i = e.nodes[0];
				if (null == i || 0 === i.terms.length) break;
				n.push(i), t = e.remainingInput
			}
			return e[i] = n
		}
	})(),
	parseExpression = (() => {
		const e = /^(\-\-|[a-z\u0240-\uffff])/i,
			t = /^([\*\+\/]|[\-]\s)/i,
			i = /^[\),]/;
		return n => {
			const r = [];
			for (; n.length && (n = n.trim(), !i.test(n));)
				if ("(" === n[0]) {
					const {
						nodes: e,
						remainingInput: t
					} = parseFunctionArguments(n);
					n = t, r.push({
						type: "function",
						name: {
							type: "ident",
							value: "calc"
						},
						arguments: e
					})
				} else if (e.test(n)) {
				const e = parseIdent(n),
					t = e.nodes[0];
				if ("(" === (n = e.remainingInput)[0]) {
					const {
						nodes: e,
						remainingInput: i
					} = parseFunctionArguments(n);
					r.push({
						type: "function",
						name: t,
						arguments: e
					}), n = i
				} else r.push(t)
			} else if (t.test(n)) r.push({
				type: "operator",
				value: n[0]
			}), n = n.slice(1);
			else {
				const {
					nodes: e,
					remainingInput: t
				} = "#" === n[0] ? parseHex(n) : parseNumber(n);
				if (0 === e.length) break;
				r.push(e[0]), n = t
			}
			return {
				nodes: [{
					type: "expression",
					terms: r
				}],
				remainingInput: n
			}
		}
	})(),
	parseIdent = (() => {
		const e = /[^a-z^0-9^_^\-^\u0240-\uffff]/i;
		return t => {
			const i = t.match(e);
			return {
				nodes: [{
					type: "ident",
					value: null == i ? t : t.substr(0, i.index)
				}],
				remainingInput: null == i ? "" : t.substr(i.index)
			}
		}
	})(),
	parseNumber = (() => {
		const e = /[\+\-]?(\d+[\.]\d+|\d+|[\.]\d+)([eE][\+\-]?\d+)?/,
			t = /^[a-z%]+/i,
			i = /^(m|mm|cm|rad|deg|[%])$/;
		return n => {
			const r = n.match(e),
				a = null == r ? "0" : r[0],
				s = (n = null == a ? n : n.slice(a.length)).match(t);
			let o = null != s && "" !== s[0] ? s[0] : null;
			const l = null == s ? n : n.slice(o.length);
			return null == o || i.test(o) || (o = null), {
				nodes: [{
					type: "number",
					number: parseFloat(a) || 0,
					unit: o
				}],
				remainingInput: l
			}
		}
	})(),
	parseHex = (() => {
		const e = /^[a-f0-9]*/i;
		return t => {
			const i = (t = t.slice(1).trim()).match(e);
			return {
				nodes: null == i ? [] : [{
					type: "hex",
					value: i[0]
				}],
				remainingInput: null == i ? t : t.slice(i[0].length)
			}
		}
	})(),
	parseFunctionArguments = e => {
		const t = [];
		for (e = e.slice(1).trim(); e.length;) {
			const i = parseExpression(e);
			if (t.push(i.nodes[0]), "," === (e = i.remainingInput.trim())[0]) e = e.slice(1).trim();
			else if (")" === e[0]) {
				e = e.slice(1);
				break
			}
		}
		return {
			nodes: t,
			remainingInput: e
		}
	},
	$visitedTypes = Symbol("visitedTypes");
class ASTWalker {
	constructor(e) {
		this[$visitedTypes] = e
	}
	walk(e, t) {
		const i = e.slice();
		for (; i.length;) {
			const e = i.shift();
			switch (this[$visitedTypes].indexOf(e.type) > -1 && t(e), e.type) {
				case "expression":
					i.unshift(...e.terms);
					break;
				case "function":
					i.unshift(e.name, ...e.arguments)
			}
		}
	}
}
const ZERO = Object.freeze({
		type: "number",
		number: 0,
		unit: null
	}),
	degreesToRadians = (e, t = 0) => {
		let {
			number: i,
			unit: n
		} = e;
		if (isFinite(i)) {
			if ("rad" === e.unit || null == e.unit) return e
		} else i = t, n = "rad";
		return {
			type: "number",
			number: ("deg" === n && null != i ? i : 0) * Math.PI / 180,
			unit: "rad"
		}
	},
	lengthToBaseMeters = (e, t = 0) => {
		let i, {
			number: n,
			unit: r
		} = e;
		if (isFinite(n)) {
			if ("m" === e.unit) return e
		} else n = t, r = "m";
		switch (r) {
			default:
				i = 1;
				break;
			case "cm":
				i = .01;
				break;
			case "mm":
				i = .001
		}
		return {
			type: "number",
			number: i * n,
			unit: "m"
		}
	},
	normalizeUnit = (() => {
		const e = e => e,
			t = {
				rad: e,
				deg: degreesToRadians,
				m: e,
				mm: lengthToBaseMeters,
				cm: lengthToBaseMeters
			};
		return (e, i = ZERO) => {
			let {
				number: n,
				unit: r
			} = e;
			if (isFinite(n) || (n = i.number, r = i.unit), null == r) return e;
			const a = t[r];
			return null == a ? i : a(e)
		}
	})();
class Hotspot extends CSS2DObject {
	constructor(e) {
		super(document.createElement("div")), this.normal = new Vector3(0, 1, 0), this.initialized = !1, this.referenceCount = 1, this.pivot = document.createElement("div"), this.slot = document.createElement("slot"), this.element.classList.add("annotation-wrapper"), this.slot.name = e.name, this.element.appendChild(this.pivot), this.pivot.appendChild(this.slot), this.updatePosition(e.position), this.updateNormal(e.normal)
	}
	get facingCamera() {
		return !this.element.classList.contains("hide")
	}
	show() {
		this.facingCamera && this.initialized || this.updateVisibility(!0)
	}
	hide() {
		!this.facingCamera && this.initialized || this.updateVisibility(!1)
	}
	increment() {
		this.referenceCount++
	}
	decrement() {
		return this.referenceCount > 0 && --this.referenceCount, 0 === this.referenceCount
	}
	updatePosition(e) {
		if (null == e) return;
		const t = parseExpressions(e)[0].terms;
		for (let e = 0; e < 3; ++e) this.position.setComponent(e, normalizeUnit(t[e]).number);
		this.updateMatrixWorld()
	}
	updateNormal(e) {
		if (null == e) return;
		const t = parseExpressions(e)[0].terms;
		for (let e = 0; e < 3; ++e) this.normal.setComponent(e, normalizeUnit(t[e]).number)
	}
	orient(e) {
		this.pivot.style.transform = `rotate(${e}rad)`
	}
	updateVisibility(e) {
		e ? this.element.classList.remove("hide") : this.element.classList.add("hide"), this.slot.assignedNodes().forEach(t => {
			if (t.nodeType !== Node.ELEMENT_NODE) return;
			const i = t,
				n = i.dataset.visibilityAttribute;
			if (null != n) {
				const t = "data-" + n;
				e ? i.setAttribute(t, "") : i.removeAttribute(t)
			}
			i.dispatchEvent(new CustomEvent("hotspot-visibility", {
				detail: {
					visible: e
				}
			}))
		}), this.initialized = !0
	}
}
const reduceVertices = (e, t, i) => {
		let n = i;
		const r = new Vector3;
		return e.traverse(e => {
			let i, a;
			e.updateWorldMatrix(!1, !1);
			const s = e.geometry;
			if (void 0 !== s)
				if (s.isGeometry) {
					const o = s.vertices;
					for (i = 0, a = o.length; i < a; i++) r.copy(o[i]), r.applyMatrix4(e.matrixWorld), n = t(n, r)
				} else if (s.isBufferGeometry) {
				const {
					position: o
				} = s.attributes;
				if (void 0 !== o)
					for (i = 0, a = o.count; i < a; i++) r.fromBufferAttribute(o, i).applyMatrix4(e.matrixWorld), n = t(n, r)
			}
		}), n
	},
	OFFSET = .002,
	LOG_MAX_RESOLUTION = 9,
	LOG_MIN_RESOLUTION = 6,
	ANIMATION_SCALING = 2;
class Shadow extends DirectionalLight {
	constructor(e, t, i) {
		super(), this.shadowMaterial = new ShadowMaterial, this.boundingBox = new Box3, this.size = new Vector3, this.isAnimated = !1, this.side = "bottom", this.needsUpdate = !1, this.intensity = 0, this.castShadow = !0, this.frustumCulled = !1, this.floor = new Mesh(new PlaneGeometry, this.shadowMaterial), this.floor.rotateX(-Math.PI / 2), this.floor.receiveShadow = !0, this.floor.castShadow = !1, this.floor.frustumCulled = !1, this.add(this.floor), e.target.add(this), this.target = e.target, this.setScene(e, t, i)
	}
	setScene(e, t, i) {
		if (this.side = i, this.isAnimated = e.animationNames.length > 0, this.boundingBox.copy(e.boundingBox), this.size.copy(e.size), "back" === this.side) {
			const {
				min: e,
				max: t
			} = this.boundingBox;
			[e.y, e.z] = [e.z, e.y], [t.y, t.z] = [t.z, t.y], [this.size.y, this.size.z] = [this.size.z, this.size.y], this.rotation.x = Math.PI / 2, this.rotation.y = Math.PI
		}
		const {
			boundingBox: n,
			size: r
		} = this;
		if (this.isAnimated) {
			const e = 2 * Math.max(r.x, r.y, r.z);
			r.y = e, n.expandByVector(r.subScalar(e).multiplyScalar(-.5)), n.max.y = n.min.y + e, r.set(e, e, e)
		}
		n.getCenter(this.floor.position);
		const a = n.max.y + .002 * r.y;
		"bottom" === i ? (this.position.y = a, this.shadow.camera.up.set(0, 0, 1)) : (this.position.y = 0, this.position.z = a, this.shadow.camera.up.set(0, 1, 0)), this.setSoftness(t)
	}
	setSoftness(e) {
		const t = Math.pow(2, 9 - 3 * e);
		this.setMapSize(t)
	}
	setMapSize(e) {
		const {
			camera: t,
			mapSize: i,
			map: n
		} = this.shadow, {
			size: r,
			boundingBox: a
		} = this;
		null != n && (n.dispose(), this.shadow.map = null), this.isAnimated && (e *= 2);
		const s = Math.floor(r.x > r.z ? e : e * r.x / r.z),
			o = Math.floor(r.x > r.z ? e * r.z / r.x : e);
		i.set(s, o);
		const l = 2.5 * r.x / s,
			c = 2.5 * r.z / o;
		t.left = -a.max.x - l, t.right = -a.min.x + l, t.bottom = a.min.z - c, t.top = a.max.z + c, this.setScaleAndOffset(t.zoom, 0), this.shadow.updateMatrices(this), this.floor.scale.set(r.x + 2 * l, r.z + 2 * c, 1), this.needsUpdate = !0
	}
	setIntensity(e) {
		this.shadowMaterial.opacity = e, e > 0 ? (this.visible = !0, this.floor.visible = !0) : (this.visible = !1, this.floor.visible = !1)
	}
	getIntensity() {
		return this.shadowMaterial.opacity
	}
	setRotation(e) {
		"bottom" === this.side ? (this.shadow.camera.up.set(Math.sin(e), 0, Math.cos(e)), this.shadow.updateMatrices(this)) : this.shadow.updateMatrices(this)
	}
	setScaleAndOffset(e, t) {
		const i = this.size.y,
			n = 1 / e,
			r = .002 * i;
		this.floor.position.y = 2 * r - i + t * n;
		const {
			camera: a
		} = this.shadow;
		a.zoom = e, a.near = 0, a.far = i * e - t, a.projectionMatrix.makeOrthographic(a.left * e, a.right * e, a.top * e, a.bottom * e, a.near, a.far), a.projectionMatrixInverse.copy(a.projectionMatrix).invert()
	}
}
const DEFAULT_FOV_DEG = 45,
	DEFAULT_HALF_FOV = 22.5 * Math.PI / 180,
	SAFE_RADIUS_RATIO = Math.sin(DEFAULT_HALF_FOV),
	DEFAULT_TAN_FOV = Math.tan(DEFAULT_HALF_FOV),
	view = new Vector3,
	target = new Vector3,
	normalWorld = new Vector3,
	raycaster = new Raycaster,
	vector3$1 = new Vector3;
class ModelScene extends Scene {
	constructor({
		canvas: e,
		element: t,
		width: i,
		height: n
	}) {
		super(), this.context = null, this.width = 1, this.height = 1, this.aspect = 1, this.isDirty = !1, this.renderCount = 0, this.camera = new PerspectiveCamera(45, 1, .1, 100), this.url = null, this.target = new Object3D, this.modelContainer = new Object3D, this.animationNames = [], this.boundingBox = new Box3, this.size = new Vector3, this.idealCameraDistance = 0, this.fieldOfViewAspect = 0, this.framedFieldOfView = 45, this.shadow = null, this.shadowIntensity = 0, this.shadowSoftness = 1, this.exposure = 1, this.canScale = !0, this.tightBounds = !1, this.goalTarget = new Vector3, this.targetDamperX = new Damper, this.targetDamperY = new Damper, this.targetDamperZ = new Damper, this._currentGLTF = null, this.cancelPendingSourceChange = null, this.animationsByName = new Map, this.currentAnimationAction = null, this.name = "ModelScene", this.element = t, this.canvas = e, this.camera = new PerspectiveCamera(45, 1, .1, 100), this.camera.name = "MainCamera", this.activeCamera = this.camera, this.add(this.target), this.setSize(i, n), this.target.name = "Target", this.modelContainer.name = "ModelContainer", this.target.add(this.modelContainer), this.mixer = new AnimationMixer(this.modelContainer)
	}
	createContext() {
		this.context = this.canvas.getContext("2d")
	}
	async setObject(e) {
		this.reset(), this.modelContainer.add(e), await this.setupScene()
	}
	async setSource(e, t) {
		if (!e || e === this.url) return void(t && t(1));
		let i;
		this.reset(), this.url = e, null != this.cancelPendingSourceChange && (this.cancelPendingSourceChange(), this.cancelPendingSourceChange = null);
		try {
			i = await new Promise(async (i, n) => {
				this.cancelPendingSourceChange = () => n();
				try {
					i(await this.element[$renderer].loader.load(e, this.element, t))
				} catch (e) {
					n(e)
				}
			})
		} catch (e) {
			if (null == e) return;
			throw e
		}
		this.reset(), this.url = e, this._currentGLTF = i, null != i && this.modelContainer.add(i.scene);
		const {
			animations: n
		} = i, r = new Map, a = [];
		for (const e of n) r.set(e.name, e), a.push(e.name);
		this.animations = n, this.animationsByName = r, this.animationNames = a, await this.setupScene()
	}
	async setupScene() {
		this.updateBoundingBox();
		let e = null;
		!0 === this.tightBounds && (await this.element.requestUpdate("cameraTarget"), e = this.getTarget()), this.updateFraming(e), this.frameModel(), this.setShadowIntensity(this.shadowIntensity), this.isDirty = !0, this.dispatchEvent({
			type: "model-load",
			url: this.url
		})
	}
	reset() {
		this.url = null;
		const e = this._currentGLTF;
		if (null != e) {
			for (const e of this.modelContainer.children) this.modelContainer.remove(e);
			e.dispose(), this._currentGLTF = null
		}
		null != this.currentAnimationAction && (this.currentAnimationAction.stop(), this.currentAnimationAction = null), this.mixer.stopAllAction(), this.mixer.uncacheRoot(this)
	}
	get currentGLTF() {
		return this._currentGLTF
	}
	setSize(e, t) {
		this.width === e && this.height === t || (this.width = Math.max(e, 1), this.height = Math.max(t, 1), this.aspect = this.width / this.height, this.frameModel(), this.isDirty = !0)
	}
	updateBoundingBox() {
		if (this.target.remove(this.modelContainer), !0 === this.tightBounds) {
			const e = (e, t) => e.expandByPoint(t);
			this.boundingBox = reduceVertices(this.modelContainer, e, new Box3)
		} else this.boundingBox.setFromObject(this.modelContainer);
		this.boundingBox.getSize(this.size), this.target.add(this.modelContainer)
	}
	updateFraming(e = null) {
		this.target.remove(this.modelContainer), null == e && (e = this.boundingBox.getCenter(new Vector3));
		const t = Math.sqrt(reduceVertices(this.modelContainer, (t, i) => Math.max(t, e.distanceToSquared(i)), 0));
		this.idealCameraDistance = t / SAFE_RADIUS_RATIO;
		this.fieldOfViewAspect = reduceVertices(this.modelContainer, (t, i) => {
			i.sub(e);
			const n = Math.sqrt(i.x * i.x + i.z * i.z);
			return Math.max(t, n / (this.idealCameraDistance - Math.abs(i.y)))
		}, 0) / DEFAULT_TAN_FOV, this.target.add(this.modelContainer)
	}
	frameModel() {
		const e = DEFAULT_TAN_FOV * Math.max(1, this.fieldOfViewAspect / this.aspect);
		this.framedFieldOfView = 2 * Math.atan(e) * 180 / Math.PI
	}
	getSize() {
		return {
			width: this.width,
			height: this.height
		}
	}
	getCamera() {
		return this.activeCamera
	}
	setCamera(e) {
		this.activeCamera = e
	}
	setTarget(e, t, i) {
		this.goalTarget.set(-e, -t, -i)
	}
	getTarget() {
		return vector3$1.copy(this.goalTarget).multiplyScalar(-1)
	}
	jumpToGoal() {
		this.updateTarget(1e4)
	}
	updateTarget(e) {
		const t = this.goalTarget,
			i = this.target.position;
		if (!t.equals(i)) {
			const n = this.idealCameraDistance;
			let {
				x: r,
				y: a,
				z: s
			} = i;
			r = this.targetDamperX.update(r, t.x, e, n), a = this.targetDamperY.update(a, t.y, e, n), s = this.targetDamperZ.update(s, t.z, e, n), this.target.position.set(r, a, s), this.target.updateMatrixWorld(), this.setShadowRotation(this.yaw), this.isDirty = !0
		}
	}
	pointTowards(e, t) {
		const {
			x: i,
			z: n
		} = this.position;
		this.yaw = Math.atan2(e - i, t - n)
	}
	set yaw(e) {
		this.rotation.y = e, this.updateMatrixWorld(!0), this.setShadowRotation(e), this.isDirty = !0
	}
	get yaw() {
		return this.rotation.y
	}
	set animationTime(e) {
		this.mixer.setTime(e)
	}
	get animationTime() {
		return null != this.currentAnimationAction ? this.currentAnimationAction.time : 0
	}
	get duration() {
		return null != this.currentAnimationAction && this.currentAnimationAction.getClip() ? this.currentAnimationAction.getClip().duration : 0
	}
	get hasActiveAnimation() {
		return null != this.currentAnimationAction
	}
	playAnimation(e = null, t = 0) {
		if (null == this._currentGLTF) return;
		const {
			animations: i
		} = this;
		if (null == i || 0 === i.length) return void console.warn("Cannot play animation (model does not have any animations)");
		let n = null;
		null != e && (n = this.animationsByName.get(e)), null == n && (n = i[0]);
		try {
			const {
				currentAnimationAction: e
			} = this;
			this.currentAnimationAction = this.mixer.clipAction(n, this).play(), this.currentAnimationAction.enabled = !0, null != e && this.currentAnimationAction !== e && this.currentAnimationAction.crossFadeFrom(e, t, !1)
		} catch (e) {
			console.error(e)
		}
	}
	stopAnimation() {
		null != this.currentAnimationAction && (this.currentAnimationAction.stop(), this.currentAnimationAction.reset(), this.currentAnimationAction = null), this.mixer.stopAllAction()
	}
	updateAnimation(e) {
		this.mixer.update(e)
	}
	updateShadow() {
		const e = this.shadow;
		if (null != e) {
			const t = "wall" === this.element.arPlacement ? "back" : "bottom";
			e.setScene(this, this.shadowSoftness, t)
		}
	}
	setShadowIntensity(e) {
		this.shadowIntensity = e;
		let t = this.shadow;
		const i = "wall" === this.element.arPlacement ? "back" : "bottom";
		null != t ? (t.setIntensity(e), t.setScene(this, this.shadowSoftness, i)) : e > 0 && (t = new Shadow(this, this.shadowSoftness, i), t.setIntensity(e), this.shadow = t)
	}
	setShadowSoftness(e) {
		this.shadowSoftness = e;
		const t = this.shadow;
		null != t && t.setSoftness(e)
	}
	setShadowRotation(e) {
		const t = this.shadow;
		null != t && t.setRotation(e)
	}
	isShadowDirty() {
		const e = this.shadow;
		if (null == e) return !1; {
			const {
				needsUpdate: t
			} = e;
			return e.needsUpdate = !1, t
		}
	}
	setShadowScaleAndOffset(e, t) {
		const i = this.shadow;
		null != i && i.setScaleAndOffset(e, t)
	}
	positionAndNormalFromPoint(e, t = this) {
		raycaster.setFromCamera(e, this.getCamera());
		const i = raycaster.intersectObject(t, !0);
		if (0 === i.length) return null;
		const n = i[0];
		return null == n.face ? null : (n.face.normal.applyNormalMatrix((new Matrix3).getNormalMatrix(n.object.matrixWorld)), {
			position: n.point,
			normal: n.face.normal
		})
	}
	addHotspot(e) {
		this.target.add(e)
	}
	removeHotspot(e) {
		this.target.remove(e)
	}
	forHotspots(e) {
		const {
			children: t
		} = this.target;
		for (let i = 0, n = t.length; i < n; i++) {
			const n = t[i];
			n instanceof Hotspot && e(n)
		}
	}
	updateHotspots(e) {
		this.forHotspots(t => {
			view.copy(e), target.setFromMatrixPosition(t.matrixWorld), view.sub(target), normalWorld.copy(t.normal).transformDirection(this.target.matrixWorld), view.dot(normalWorld) < 0 ? t.hide() : t.show()
		})
	}
	orientHotspots(e) {
		this.forHotspots(t => {
			t.orient(e)
		})
	}
	setHotspotsVisibility(e) {
		this.forHotspots(t => {
			t.visible = e
		})
	}
}
var _mipmapMaterial = _getMipmapMaterial(),
	_mesh = new Mesh(new PlaneGeometry(2, 2), _mipmapMaterial),
	_flatCamera = new OrthographicCamera(0, 1, 0, 1, 0, 1),
	_tempTarget = null,
	_renderer = null;

function RoughnessMipmapper(e) {
	(_renderer = e).compile(_mesh, _flatCamera)
}

function _getMipmapMaterial() {
	var e = new RawShaderMaterial({
		uniforms: {
			roughnessMap: {
				value: null
			},
			normalMap: {
				value: null
			},
			texelSize: {
				value: new Vector2(1, 1)
			}
		},
		vertexShader: "\n\t\t\tprecision mediump float;\n\t\t\tprecision mediump int;\n\n\t\t\tattribute vec3 position;\n\t\t\tattribute vec2 uv;\n\n\t\t\tvarying vec2 vUv;\n\n\t\t\tvoid main() {\n\n\t\t\t\tvUv = uv;\n\n\t\t\t\tgl_Position = vec4( position, 1.0 );\n\n\t\t\t}\n\t\t",
		fragmentShader: "\n\t\t\tprecision mediump float;\n\t\t\tprecision mediump int;\n\n\t\t\tvarying vec2 vUv;\n\n\t\t\tuniform sampler2D roughnessMap;\n\t\t\tuniform sampler2D normalMap;\n\t\t\tuniform vec2 texelSize;\n\n\t\t\t#define ENVMAP_TYPE_CUBE_UV\n\n\t\t\tvec4 envMapTexelToLinear( vec4 a ) { return a; }\n\n\t\t\t#include <cube_uv_reflection_fragment>\n\n\t\t\tfloat roughnessToVariance( float roughness ) {\n\n\t\t\t\tfloat variance = 0.0;\n\n\t\t\t\tif ( roughness >= r1 ) {\n\n\t\t\t\t\tvariance = ( r0 - roughness ) * ( v1 - v0 ) / ( r0 - r1 ) + v0;\n\n\t\t\t\t} else if ( roughness >= r4 ) {\n\n\t\t\t\t\tvariance = ( r1 - roughness ) * ( v4 - v1 ) / ( r1 - r4 ) + v1;\n\n\t\t\t\t} else if ( roughness >= r5 ) {\n\n\t\t\t\t\tvariance = ( r4 - roughness ) * ( v5 - v4 ) / ( r4 - r5 ) + v4;\n\n\t\t\t\t} else {\n\n\t\t\t\t\tfloat roughness2 = roughness * roughness;\n\n\t\t\t\t\tvariance = 1.79 * roughness2 * roughness2;\n\n\t\t\t\t}\n\n\t\t\t\treturn variance;\n\n\t\t\t}\n\n\t\t\tfloat varianceToRoughness( float variance ) {\n\n\t\t\t\tfloat roughness = 0.0;\n\n\t\t\t\tif ( variance >= v1 ) {\n\n\t\t\t\t\troughness = ( v0 - variance ) * ( r1 - r0 ) / ( v0 - v1 ) + r0;\n\n\t\t\t\t} else if ( variance >= v4 ) {\n\n\t\t\t\t\troughness = ( v1 - variance ) * ( r4 - r1 ) / ( v1 - v4 ) + r1;\n\n\t\t\t\t} else if ( variance >= v5 ) {\n\n\t\t\t\t\troughness = ( v4 - variance ) * ( r5 - r4 ) / ( v4 - v5 ) + r4;\n\n\t\t\t\t} else {\n\n\t\t\t\t\troughness = pow( 0.559 * variance, 0.25 ); // 0.559 = 1.0 / 1.79\n\n\t\t\t\t}\n\n\t\t\t\treturn roughness;\n\n\t\t\t}\n\n\t\t\tvoid main() {\n\n\t\t\t\tgl_FragColor = texture2D( roughnessMap, vUv, - 1.0 );\n\n\t\t\t\tif ( texelSize.x == 0.0 ) return;\n\n\t\t\t\tfloat roughness = gl_FragColor.g;\n\n\t\t\t\tfloat variance = roughnessToVariance( roughness );\n\n\t\t\t\tvec3 avgNormal;\n\n\t\t\t\tfor ( float x = - 1.0; x < 2.0; x += 2.0 ) {\n\n\t\t\t\t\tfor ( float y = - 1.0; y < 2.0; y += 2.0 ) {\n\n\t\t\t\t\t\tvec2 uv = vUv + vec2( x, y ) * 0.25 * texelSize;\n\n\t\t\t\t\t\tavgNormal += normalize( texture2D( normalMap, uv, - 1.0 ).xyz - 0.5 );\n\n\t\t\t\t\t}\n\n\t\t\t\t}\n\n\t\t\t\tvariance += 1.0 - 0.25 * length( avgNormal );\n\n\t\t\t\tgl_FragColor.g = varianceToRoughness( variance );\n\n\t\t\t}\n\t\t",
		blending: 0,
		depthTest: !1,
		depthWrite: !1
	});
	return e.type = "RoughnessMipmapper", e
}
RoughnessMipmapper.prototype = {
	constructor: RoughnessMipmapper,
	generateMipmaps: function(e) {
		if ("roughnessMap" in e != !1) {
			var {
				roughnessMap: t,
				normalMap: i
			} = e;
			if (null !== t && null !== i && t.generateMipmaps && !e.userData.roughnessUpdated) {
				e.userData.roughnessUpdated = !0;
				var n = Math.max(t.image.width, i.image.width),
					r = Math.max(t.image.height, i.image.height);
				if (MathUtils.isPowerOfTwo(n) && MathUtils.isPowerOfTwo(r)) {
					var a = _renderer.getRenderTarget(),
						s = _renderer.autoClear;
					if (_renderer.autoClear = !1, null !== _tempTarget && _tempTarget.width === n && _tempTarget.height === r || (null !== _tempTarget && _tempTarget.dispose(), (_tempTarget = new WebGLRenderTarget(n, r, {
							depthBuffer: !1
						})).scissorTest = !0), n !== t.image.width || r !== t.image.height) {
						var o = {
								wrapS: t.wrapS,
								wrapT: t.wrapT,
								magFilter: t.magFilter,
								minFilter: t.minFilter,
								depthBuffer: !1
							},
							l = new WebGLRenderTarget(n, r, o);
						l.texture.generateMipmaps = !0, _renderer.setRenderTarget(l), e.roughnessMap = l.texture, e.metalnessMap == t && (e.metalnessMap = e.roughnessMap), e.aoMap == t && (e.aoMap = e.roughnessMap)
					}
					_mipmapMaterial.uniforms.roughnessMap.value = t, _mipmapMaterial.uniforms.normalMap.value = i;
					for (var c = new Vector2(0, 0), h = _mipmapMaterial.uniforms.texelSize.value, u = 0; n >= 1 && r >= 1; ++u, n /= 2, r /= 2) h.set(1 / n, 1 / r), 0 == u && h.set(0, 0), _tempTarget.viewport.set(c.x, c.y, n, r), _tempTarget.scissor.set(c.x, c.y, n, r), _renderer.setRenderTarget(_tempTarget), _renderer.render(_mesh, _flatCamera), _renderer.copyFramebufferToTexture(c, e.roughnessMap, u), _mipmapMaterial.uniforms.roughnessMap.value = e.roughnessMap;
					t !== e.roughnessMap && t.dispose(), _renderer.setRenderTarget(a), _renderer.autoClear = s
				}
			}
		}
	},
	dispose: function() {
		_mipmapMaterial.dispose(), _mesh.geometry.dispose(), null != _tempTarget && _tempTarget.dispose()
	}
};
const deserializeUrl = e => e && "null" !== e ? toFullUrl(e) : null,
	assertIsArCandidate = () => {
		if (IS_WEBXR_AR_CANDIDATE) return;
		const e = [];
		throw HAS_WEBXR_DEVICE_API || e.push("WebXR Device API"), HAS_WEBXR_HIT_TEST_API || e.push("WebXR Hit Test API"), new Error("The following APIs are required for AR, but are missing in this browser: " + e.join(", "))
	},
	toFullUrl = e => new URL(e, window.location.toString()).toString(),
	throttle = (e, t) => {
		let i = null;
		const n = (...n) => {
			null == i && (e(...n), i = self.setTimeout(() => i = null, t))
		};
		return n.flush = () => {
			null != i && (self.clearTimeout(i), i = null)
		}, n
	},
	debounce = (e, t) => {
		let i = null;
		return (...n) => {
			null != i && self.clearTimeout(i), i = self.setTimeout(() => {
				i = null, e(...n)
			}, t)
		}
	},
	clamp = (e, t, i) => Math.max(t, Math.min(i, e)),
	CAPPED_DEVICE_PIXEL_RATIO = 1,
	resolveDpr = (() => {
		const e = (() => {
			const e = null != document.head ? Array.from(document.head.querySelectorAll("meta")) : [];
			for (const t of e)
				if ("viewport" === t.name) return !0;
			return !1
		})();
		return e || console.warn('No <meta name="viewport"> detected; <model-viewer> will cap pixel density at 1.'), () => e ? window.devicePixelRatio : 1
	})(),
	isDebugMode = (() => {
		const e = new RegExp("[?&]model-viewer-debug-mode(&|$)");
		return () => self.ModelViewerElement && self.ModelViewerElement.debugMode || self.location && self.location.search && self.location.search.match(e)
	})(),
	getFirstMapKey = e => {
		if (null != e.keys) return e.keys().next().value || null;
		let t = null;
		try {
			e.forEach((e, i, n) => {
				throw t = i, new Error
			})
		} catch (e) {}
		return t
	},
	timePasses = (e = 0) => new Promise(t => setTimeout(t, e)),
	waitForEvent = (e, t, i = null) => new Promise(n => {
		e.addEventListener(t, (function r(a) {
			i && !i(a) || (n(a), e.removeEventListener(t, r))
		}))
	}),
	RADIUS = .2,
	LINE_WIDTH = .03,
	MAX_OPACITY = .75,
	SEGMENTS = 12,
	DELTA_PHI = Math.PI / 24,
	vector2 = new Vector2,
	addCorner = (e, t, i) => {
		let n = t > 0 ? i > 0 ? 0 : -Math.PI / 2 : i > 0 ? Math.PI / 2 : Math.PI;
		for (let r = 0; r <= 12; ++r) e.push(t + .17 * Math.cos(n), i + .17 * Math.sin(n), 0, t + .2 * Math.cos(n), i + .2 * Math.sin(n), 0), n += DELTA_PHI
	};
class PlacementBox extends Mesh {
	constructor(e, t) {
		const i = new BufferGeometry,
			n = [],
			r = [],
			{
				size: a,
				boundingBox: s
			} = e,
			o = a.x / 2,
			l = ("back" === t ? a.y : a.z) / 2;
		addCorner(r, o, l), addCorner(r, -o, l), addCorner(r, -o, -l), addCorner(r, o, -l);
		const c = r.length / 3;
		for (let e = 0; e < c - 2; e += 2) n.push(e, e + 1, e + 3, e, e + 3, e + 2);
		const h = c - 2;
		n.push(h, h + 1, 1, h, 1, 0), i.setAttribute("position", new Float32BufferAttribute(r, 3)), i.setIndex(n), super(i), this.side = t;
		const u = this.material;
		switch (u.side = 2, u.transparent = !0, u.opacity = 0, this.goalOpacity = 0, this.opacityDamper = new Damper, this.hitPlane = new Mesh(new PlaneGeometry(2 * (o + .2), 2 * (l + .2))), this.hitPlane.visible = !1, this.add(this.hitPlane), s.getCenter(this.position), t) {
			case "bottom":
				this.rotateX(-Math.PI / 2), this.shadowHeight = s.min.y, this.position.y = this.shadowHeight;
				break;
			case "back":
				this.shadowHeight = s.min.z, this.position.z = this.shadowHeight
		}
		e.target.add(this)
	}
	getHit(e, t, i) {
		vector2.set(t, -i), this.hitPlane.visible = !0;
		const n = e.positionAndNormalFromPoint(vector2, this.hitPlane);
		return this.hitPlane.visible = !1, null == n ? null : n.position
	}
	set offsetHeight(e) {
		"back" === this.side ? this.position.z = this.shadowHeight + e : this.position.y = this.shadowHeight + e
	}
	get offsetHeight() {
		return "back" === this.side ? this.position.z - this.shadowHeight : this.position.y - this.shadowHeight
	}
	set show(e) {
		this.goalOpacity = e ? .75 : 0
	}
	updateOpacity(e) {
		const t = this.material;
		t.opacity = this.opacityDamper.update(t.opacity, this.goalOpacity, e, 1), this.visible = t.opacity > 0
	}
	dispose() {
		const {
			geometry: e,
			material: t
		} = this.hitPlane;
		e.dispose(), t.dispose(), this.geometry.dispose(), this.material.dispose()
	}
}
const AR_SHADOW_INTENSITY = .3,
	ROTATION_RATE = 1.5,
	HIT_ANGLE_DEG = 20,
	INTRO_DAMPER_RATE = .4,
	SCALE_SNAP_HIGH = 1.2,
	SCALE_SNAP_LOW = 1 / 1.2,
	MIN_VIEWPORT_SCALE = .25,
	ARStatus = {
		NOT_PRESENTING: "not-presenting",
		SESSION_STARTED: "session-started",
		OBJECT_PLACED: "object-placed",
		FAILED: "failed"
	},
	vector3 = new Vector3,
	matrix4 = new Matrix4,
	hitPosition = new Vector3;
class ARRenderer extends EventDispatcher {
	constructor(e) {
		super(), this.renderer = e, this.camera = new PerspectiveCamera, this.currentSession = null, this.placeOnWall = !1, this.placementBox = null, this.lastTick = null, this.turntableRotation = null, this.oldShadowIntensity = null, this.oldBackground = null, this.rafId = null, this.refSpace = null, this.viewerRefSpace = null, this.frame = null, this.initialHitSource = null, this.transientHitTestSource = null, this.inputSource = null, this._presentedScene = null, this.resolveCleanup = null, this.exitWebXRButtonContainer = null, this.initialModelToWorld = null, this.initialized = !1, this.oldTarget = new Vector3, this.placementComplete = !1, this.isTranslating = !1, this.isRotating = !1, this.isScaling = !1, this.lastDragPosition = new Vector3, this.lastScalar = 0, this.goalPosition = new Vector3, this.goalYaw = 0, this.goalScale = 1, this.xDamper = new Damper, this.yDamper = new Damper, this.zDamper = new Damper, this.yawDamper = new Damper, this.scaleDamper = new Damper, this.damperRate = 1, this.onExitWebXRButtonContainerClick = () => this.stopPresenting(), this.onUpdateScene = () => {
			null != this.placementBox && this.isPresenting && (this.placementBox.dispose(), this.placementBox = new PlacementBox(this.presentedScene, this.placeOnWall ? "back" : "bottom"))
		}, this.onSelectStart = e => {
			const t = this.transientHitTestSource;
			if (null == t) return;
			const i = this.frame.getHitTestResultsForTransientInput(t),
				n = this.presentedScene,
				r = this.placementBox;
			if (1 === i.length) {
				this.inputSource = e.inputSource;
				const {
					axes: t
				} = this.inputSource.gamepad, i = r.getHit(this.presentedScene, t[0], t[1]);
				r.show = !0, null != i ? (this.isTranslating = !0, this.lastDragPosition.copy(i)) : !1 === this.placeOnWall && (this.isRotating = !0, this.lastScalar = t[0])
			} else 2 === i.length && n.canScale && (r.show = !0, this.isScaling = !0, this.lastScalar = this.fingerSeparation(i) / n.scale.x)
		}, this.onSelectEnd = () => {
			this.isTranslating = !1, this.isRotating = !1, this.isScaling = !1, this.inputSource = null, this.goalPosition.y += this.placementBox.offsetHeight * this.presentedScene.scale.x, this.placementBox.show = !1
		}, this.threeRenderer = e.threeRenderer, this.camera.matrixAutoUpdate = !1
	}
	async resolveARSession(e) {
		assertIsArCandidate();
		const t = await navigator.xr.requestSession("immersive-ar", {
				requiredFeatures: ["hit-test"],
				optionalFeatures: ["dom-overlay"],
				domOverlay: {
					root: e.element.shadowRoot.querySelector("div.default")
				}
			}),
			i = this.threeRenderer.getContext();
		await i.makeXRCompatible(), t.updateRenderState({
			baseLayer: new XRWebGLLayer(t, i, {
				alpha: !0
			})
		});
		let n = new Promise((e, i) => {
			t.requestAnimationFrame(() => e())
		});
		await n, e.element[$onResize](window.screen);
		const {
			framebuffer: r,
			framebufferWidth: a,
			framebufferHeight: s
		} = t.renderState.baseLayer;
		this.threeRenderer.setFramebuffer(r), this.threeRenderer.setPixelRatio(1), this.threeRenderer.setSize(a, s, !1);
		const o = e.element.shadowRoot.querySelector(".slot.exit-webxr-ar-button");
		return o.classList.add("enabled"), o.addEventListener("click", this.onExitWebXRButtonContainerClick), this.exitWebXRButtonContainer = o, t
	}
	get presentedScene() {
		return this._presentedScene
	}
	async supportsPresentation() {
		try {
			return assertIsArCandidate(), await navigator.xr.isSessionSupported("immersive-ar")
		} catch (e) {
			return !1
		}
	}
	async present(e) {
		this.isPresenting && console.warn("Cannot present while a model is already presenting");
		let t = new Promise((e, t) => {
			requestAnimationFrame(() => e())
		});
		e.setHotspotsVisibility(!1), e.isDirty = !0, await t, this._presentedScene = e;
		const i = await this.resolveARSession(e);
		i.addEventListener("end", () => {
			this.postSessionCleanup()
		}, {
			once: !0
		}), this.refSpace = await i.requestReferenceSpace("local"), this.viewerRefSpace = await i.requestReferenceSpace("viewer"), e.setCamera(this.camera), this.initialized = !1, this.damperRate = .4, this.turntableRotation = e.yaw, e.yaw = 0, this.goalYaw = 0, this.goalScale = 1, this.oldBackground = e.background, e.background = null, this.oldShadowIntensity = e.shadowIntensity, e.setShadowIntensity(0), this.oldTarget.copy(e.getTarget()), e.addEventListener("model-load", this.onUpdateScene);
		const n = 20 * Math.PI / 180,
			r = !0 === this.placeOnWall ? void 0 : new XRRay(new DOMPoint(0, 0, 0), {
				x: 0,
				y: -Math.sin(n),
				z: -Math.cos(n)
			});
		i.requestHitTestSource({
			space: this.viewerRefSpace,
			offsetRay: r
		}).then(e => {
			this.initialHitSource = e
		}), this.currentSession = i, this.placementBox = new PlacementBox(e, this.placeOnWall ? "back" : "bottom"), this.placementComplete = !1, this.lastTick = performance.now(), this.tick()
	}
	async stopPresenting() {
		if (!this.isPresenting) return;
		const e = new Promise(e => {
			this.resolveCleanup = e
		});
		try {
			await this.currentSession.end(), await e
		} catch (e) {
			console.warn("Error while trying to end AR session"), console.warn(e), this.postSessionCleanup()
		}
	}
	get isPresenting() {
		return null != this.presentedScene
	}
	get target() {
		return this.oldTarget
	}
	updateTarget() {
		const e = this.presentedScene;
		if (null != e) {
			const t = e.getTarget();
			this.oldTarget.copy(t), this.placeOnWall ? e.setTarget(t.x, t.y, e.boundingBox.min.z) : e.setTarget(t.x, e.boundingBox.min.y, t.z)
		}
	}
	postSessionCleanup() {
		this.threeRenderer.setFramebuffer(null);
		const e = this.currentSession;
		null != e && (e.removeEventListener("selectstart", this.onSelectStart), e.removeEventListener("selectend", this.onSelectEnd), e.cancelAnimationFrame(this.rafId), this.currentSession = null);
		const t = this.presentedScene;
		if (null != t) {
			const {
				target: e,
				element: i
			} = t;
			t.setCamera(t.camera), e.remove(this.placementBox), t.position.set(0, 0, 0), t.scale.set(1, 1, 1), t.setShadowScaleAndOffset(1, 0);
			const n = this.turntableRotation;
			null != n && (t.yaw = n);
			const r = this.oldShadowIntensity;
			null != r && t.setShadowIntensity(r);
			const a = this.oldBackground;
			null != a && (t.background = a);
			const s = this.oldTarget;
			t.setTarget(s.x, s.y, s.z), t.removeEventListener("model-load", this.onUpdateScene), t.orientHotspots(0), i.requestUpdate("cameraTarget"), i[$onResize](i.getBoundingClientRect())
		}
		this.renderer.height = 0;
		const i = this.exitWebXRButtonContainer;
		null != i && (i.classList.remove("enabled"), i.removeEventListener("click", this.onExitWebXRButtonContainerClick), this.exitWebXRButtonContainer = null);
		const n = this.transientHitTestSource;
		null != n && (n.cancel(), this.transientHitTestSource = null);
		const r = this.initialHitSource;
		null != r && (r.cancel(), this.initialHitSource = null), null != this.placementBox && (this.placementBox.dispose(), this.placementBox = null), this.lastTick = null, this.turntableRotation = null, this.oldShadowIntensity = null, this.oldBackground = null, this.rafId = null, this.refSpace = null, this._presentedScene = null, this.viewerRefSpace = null, this.frame = null, this.inputSource = null, null != this.resolveCleanup && this.resolveCleanup(), this.dispatchEvent({
			type: "status",
			status: ARStatus.NOT_PRESENTING
		})
	}
	updateCamera(e) {
		const {
			camera: t
		} = this, {
			matrix: i
		} = t;
		if (i.fromArray(e.transform.matrix), t.updateMatrixWorld(!0), t.position.setFromMatrixPosition(i), !this.initialized) {
			t.projectionMatrix.fromArray(e.projectionMatrix), t.projectionMatrixInverse.copy(t.projectionMatrix).invert();
			const i = this.presentedScene;
			t.getWorldDirection(vector3), i.yaw = Math.atan2(-vector3.x, -vector3.z), this.goalYaw = i.yaw, this.initialModelToWorld = (new Matrix4).copy(i.matrixWorld), i.setHotspotsVisibility(!0), this.initialized = !0, this.dispatchEvent({
				type: "status",
				status: ARStatus.SESSION_STARTED
			})
		}
		if (null != this.initialHitSource) {
			const {
				position: e,
				idealCameraDistance: i
			} = this.presentedScene;
			t.getWorldDirection(e), e.multiplyScalar(i), e.add(t.position)
		}
		if (e.requestViewportScale && e.recommendedViewportScale) {
			const t = e.recommendedViewportScale;
			e.requestViewportScale(Math.max(t, .25))
		}
		const n = this.currentSession.renderState.baseLayer.getViewport(e);
		this.threeRenderer.setViewport(n.x, n.y, n.width, n.height), this.presentedScene.orientHotspots(Math.atan2(i.elements[1], i.elements[5]))
	}
	placeInitially(e) {
		const t = this.initialHitSource;
		if (null == t) return;
		const i = e.getHitTestResults(t);
		if (0 == i.length) return;
		const n = i[0],
			r = this.getHitPoint(n);
		if (null == r) return;
		this.placeModel(r), t.cancel(), this.initialHitSource = null;
		const {
			session: a
		} = e;
		a.addEventListener("selectstart", this.onSelectStart), a.addEventListener("selectend", this.onSelectEnd), a.requestHitTestSourceForTransientInput({
			profile: "generic-touchscreen"
		}).then(e => {
			this.transientHitTestSource = e
		})
	}
	getHitPoint(e) {
		const t = e.getPose(this.refSpace);
		if (null == t) return null;
		const i = matrix4.fromArray(t.transform.matrix);
		return !0 === this.placeOnWall && (this.goalYaw = Math.atan2(i.elements[4], i.elements[6])), i.elements[5] > .75 !== this.placeOnWall ? hitPosition.setFromMatrixPosition(i) : null
	}
	placeModel(e) {
		const t = this.presentedScene;
		this.placementBox.show = !0;
		const i = this.goalPosition;
		if (i.copy(e), !1 === this.placeOnWall) {
			const n = e.y,
				r = this.camera.position.clone(),
				a = e.clone().sub(r).normalize();
			r.sub(a.multiplyScalar(t.idealCameraDistance));
			const s = new Ray(r, a.normalize()),
				o = this.initialModelToWorld,
				l = (new Vector3).setFromMatrixPosition(o).add(e);
			o.setPosition(l), s.applyMatrix4(o.invert());
			const {
				max: c
			} = t.boundingBox;
			c.y += 10, s.intersectBox(t.boundingBox, l), c.y -= 10, null != l && (l.applyMatrix4(o), i.add(e).sub(l)), i.y = n
		}
		this.updateTarget(), this.dispatchEvent({
			type: "status",
			status: ARStatus.OBJECT_PLACED
		})
	}
	fingerSeparation(e) {
		const t = e[0].inputSource.gamepad.axes,
			i = e[1].inputSource.gamepad.axes,
			n = i[0] - t[0],
			r = i[1] - t[1];
		return Math.sqrt(n * n + r * r)
	}
	processInput(e) {
		const t = this.transientHitTestSource;
		if (null == t) return;
		if (!this.isTranslating && !this.isScaling && !this.isRotating) return;
		const i = e.getHitTestResultsForTransientInput(t),
			n = this.presentedScene,
			r = n.scale.x;
		if (this.isScaling)
			if (i.length < 2) this.isScaling = !1;
			else {
				const e = this.fingerSeparation(i) / this.lastScalar;
				this.goalScale = e < 1.2 && e > 1 / 1.2 ? 1 : e
			}
		else {
			if (2 === i.length && n.canScale) return this.isTranslating = !1, this.isRotating = !1, this.isScaling = !0, void(this.lastScalar = this.fingerSeparation(i) / r);
			if (this.isRotating) {
				const e = this.inputSource.gamepad.axes[0];
				this.goalYaw += 1.5 * (e - this.lastScalar), this.lastScalar = e
			} else this.isTranslating && i.forEach(e => {
				if (e.inputSource !== this.inputSource || e.results.length < 1) return;
				const t = this.getHitPoint(e.results[0]);
				if (null != t) {
					if (this.goalPosition.sub(this.lastDragPosition), !1 === this.placeOnWall) {
						const e = t.y - this.lastDragPosition.y;
						if (e < 0) {
							this.placementBox.offsetHeight = e / r, this.presentedScene.setShadowScaleAndOffset(r, e);
							const i = vector3.copy(this.camera.position),
								n = -e / (i.y - t.y);
							i.multiplyScalar(n), t.multiplyScalar(1 - n).add(i)
						}
					}
					this.goalPosition.add(t), this.lastDragPosition.copy(t)
				}
			})
		}
	}
	moveScene(e) {
		const t = this.presentedScene,
			{
				position: i,
				yaw: n,
				idealCameraDistance: r
			} = t,
			a = this.goalPosition,
			s = t.scale.x,
			o = this.placementBox;
		if (null == this.initialHitSource && (!a.equals(i) || this.goalScale !== s)) {
			let {
				x: n,
				y: l,
				z: c
			} = i;
			e *= this.damperRate, n = this.xDamper.update(n, a.x, e, r), l = this.yDamper.update(l, a.y, e, r), c = this.zDamper.update(c, a.z, e, r), i.set(n, l, c);
			const h = this.scaleDamper.update(s, this.goalScale, e, 1);
			if (t.scale.set(h, h, h), !this.isTranslating) {
				const e = a.y - l;
				this.placementComplete && !1 === this.placeOnWall ? (o.offsetHeight = e / h, t.setShadowScaleAndOffset(h, e)) : 0 === e && (this.placementComplete = !0, o.show = !1, t.setShadowIntensity(.3), this.damperRate = 1)
			}
		}
		o.updateOpacity(e), t.updateTarget(e), t.yaw = this.yawDamper.update(n, this.goalYaw, e, Math.PI)
	}
	tick() {
		this.rafId = this.currentSession.requestAnimationFrame((e, t) => this.onWebXRFrame(e, t))
	}
	onWebXRFrame(e, t) {
		this.frame = t;
		const i = t.getViewerPose(this.refSpace);
		this.tick();
		const n = this.presentedScene;
		if (null == i || null == n) return;
		let r = !0;
		for (const a of i.views) {
			if (this.updateCamera(a), r) {
				this.placeInitially(t), this.processInput(t);
				const i = e - this.lastTick;
				this.moveScene(i), this.renderer.preRender(n, e, i), this.lastTick = e
			}
			const i = this.threeRenderer.context;
			i.depthMask(!1), i.clear(i.DEPTH_BUFFER_BIT), i.depthMask(!0), this.threeRenderer.render(n, this.camera), r = !1
		}
	}
}
class Debugger {
	constructor(e) {
		e.threeRenderer.debug = {
			checkShaderErrors: !0
		}, Promise.resolve().then(() => {
			self.dispatchEvent(new CustomEvent("model-viewer-renderer-debug", {
				detail: {
					renderer: e,
					THREE: {
						ShaderMaterial: ShaderMaterial,
						Texture: Texture$1,
						Mesh: Mesh,
						Scene: Scene,
						PlaneBufferGeometry: PlaneGeometry,
						OrthographicCamera: OrthographicCamera,
						WebGLRenderTarget: WebGLRenderTarget
					}
				}
			}))
		})
	}
	addScene(e) {
		self.dispatchEvent(new CustomEvent("model-viewer-scene-added-debug", {
			detail: {
				scene: e
			}
		}))
	}
	removeScene(e) {
		self.dispatchEvent(new CustomEvent("model-viewer-scene-removed-debug", {
			detail: {
				scene: e
			}
		}))
	}
}
var SkeletonUtils = {
	retarget: function() {
		var e = new Vector3,
			t = new Quaternion,
			i = new Vector3,
			n = new Matrix4,
			r = new Matrix4,
			a = new Matrix4;
		return function(s, o, l) {
			(l = l || {}).preserveMatrix = void 0 === l.preserveMatrix || l.preserveMatrix, l.preservePosition = void 0 === l.preservePosition || l.preservePosition, l.preserveHipPosition = void 0 !== l.preserveHipPosition && l.preserveHipPosition, l.useTargetMatrix = void 0 !== l.useTargetMatrix && l.useTargetMatrix, l.hip = void 0 !== l.hip ? l.hip : "hip", l.names = l.names || {};
			var c, h, u, d, p, m, A = o.isObject3D ? o.skeleton.bones : this.getBones(o),
				g = s.isObject3D ? s.skeleton.bones : this.getBones(s);
			if (s.isObject3D ? s.skeleton.pose() : (l.useTargetMatrix = !0, l.preserveMatrix = !1), l.preservePosition)
				for (p = [], m = 0; m < g.length; m++) p.push(g[m].position.clone());
			if (l.preserveMatrix)
				for (s.updateMatrixWorld(), s.matrixWorld.identity(), m = 0; m < s.children.length; ++m) s.children[m].updateMatrixWorld(!0);
			if (l.offsets)
				for (c = [], m = 0; m < g.length; ++m) h = g[m], u = l.names[h.name] || h.name, l.offsets && l.offsets[u] && (h.matrix.multiply(l.offsets[u]), h.matrix.decompose(h.position, h.quaternion, h.scale), h.updateMatrixWorld()), c.push(h.matrixWorld.clone());
			for (m = 0; m < g.length; ++m) {
				if (h = g[m], u = l.names[h.name] || h.name, d = this.getBoneByName(u, A), a.copy(h.matrixWorld), d) {
					if (d.updateMatrixWorld(), l.useTargetMatrix ? r.copy(d.matrixWorld) : (r.copy(s.matrixWorld).invert(), r.multiply(d.matrixWorld)), i.setFromMatrixScale(r), r.scale(i.set(1 / i.x, 1 / i.y, 1 / i.z)), a.makeRotationFromQuaternion(t.setFromRotationMatrix(r)), s.isObject3D) {
						var f = g.indexOf(h),
							v = c ? c[f] : n.copy(s.skeleton.boneInverses[f]).invert();
						a.multiply(v)
					}
					a.copyPosition(r)
				}
				h.parent && h.parent.isBone ? (h.matrix.copy(h.parent.matrixWorld).invert(), h.matrix.multiply(a)) : h.matrix.copy(a), l.preserveHipPosition && u === l.hip && h.matrix.setPosition(e.set(0, h.position.y, 0)), h.matrix.decompose(h.position, h.quaternion, h.scale), h.updateMatrixWorld()
			}
			if (l.preservePosition)
				for (m = 0; m < g.length; ++m) h = g[m], (u = l.names[h.name] || h.name) !== l.hip && h.position.copy(p[m]);
			l.preserveMatrix && s.updateMatrixWorld(!0)
		}
	}(),
	retargetClip: function(e, t, i, n) {
		(n = n || {}).useFirstFramePosition = void 0 !== n.useFirstFramePosition && n.useFirstFramePosition, n.fps = void 0 !== n.fps ? n.fps : 30, n.names = n.names || [], t.isObject3D || (t = this.getHelperFromSkeleton(t));
		var r, a, s, o, l, c, h = Math.round(i.duration * (n.fps / 1e3) * 1e3),
			u = 1 / n.fps,
			d = [],
			p = new AnimationMixer(t),
			m = this.getBones(e.skeleton),
			A = [];
		for (p.clipAction(i).play(), p.update(0), t.updateMatrixWorld(), l = 0; l < h; ++l) {
			var g = l * u;
			for (this.retarget(e, t, n), c = 0; c < m.length; ++c) o = n.names[m[c].name] || m[c].name, this.getBoneByName(o, t.skeleton) && (a = m[c], s = A[c] = A[c] || {
				bone: a
			}, n.hip === o && (s.pos || (s.pos = {
				times: new Float32Array(h),
				values: new Float32Array(3 * h)
			}), n.useFirstFramePosition && (0 === l && (r = a.position.clone()), a.position.sub(r)), s.pos.times[l] = g, a.position.toArray(s.pos.values, 3 * l)), s.quat || (s.quat = {
				times: new Float32Array(h),
				values: new Float32Array(4 * h)
			}), s.quat.times[l] = g, a.quaternion.toArray(s.quat.values, 4 * l));
			p.update(u), t.updateMatrixWorld()
		}
		for (l = 0; l < A.length; ++l)(s = A[l]) && (s.pos && d.push(new VectorKeyframeTrack(".bones[" + s.bone.name + "].position", s.pos.times, s.pos.values)), d.push(new QuaternionKeyframeTrack(".bones[" + s.bone.name + "].quaternion", s.quat.times, s.quat.values)));
		return p.uncacheAction(i), new AnimationClip(i.name, -1, d)
	},
	getHelperFromSkeleton: function(e) {
		var t = new SkeletonHelper(e.bones[0]);
		return t.skeleton = e, t
	},
	getSkeletonOffsets: function() {
		var e = new Vector3,
			t = new Vector3,
			i = new Vector3,
			n = new Vector3,
			r = new Vector2,
			a = new Vector2;
		return function(s, o, l) {
			(l = l || {}).hip = void 0 !== l.hip ? l.hip : "hip", l.names = l.names || {}, o.isObject3D || (o = this.getHelperFromSkeleton(o));
			var c, h, u, d, p = Object.keys(l.names),
				m = Object.values(l.names),
				A = o.isObject3D ? o.skeleton.bones : this.getBones(o),
				g = s.isObject3D ? s.skeleton.bones : this.getBones(s),
				f = [];
			for (s.skeleton.pose(), d = 0; d < g.length; ++d)
				if (c = g[d], u = l.names[c.name] || c.name, (h = this.getBoneByName(u, A)) && u !== l.hip) {
					var v = this.getNearestBone(c.parent, p),
						y = this.getNearestBone(h.parent, m);
					v.updateMatrixWorld(), y.updateMatrixWorld(), e.setFromMatrixPosition(v.matrixWorld), t.setFromMatrixPosition(c.matrixWorld), i.setFromMatrixPosition(y.matrixWorld), n.setFromMatrixPosition(h.matrixWorld), r.subVectors(new Vector2(t.x, t.y), new Vector2(e.x, e.y)).normalize(), a.subVectors(new Vector2(n.x, n.y), new Vector2(i.x, i.y)).normalize();
					var E = r.angle() - a.angle(),
						_ = (new Matrix4).makeRotationFromEuler(new Euler(0, 0, E));
					c.matrix.multiply(_), c.matrix.decompose(c.position, c.quaternion, c.scale), c.updateMatrixWorld(), f[u] = _
				} return f
		}
	}(),
	renameBones: function(e, t) {
		for (var i = this.getBones(e), n = 0; n < i.length; ++n) {
			var r = i[n];
			t[r.name] && (r.name = t[r.name])
		}
		return this
	},
	getBones: function(e) {
		return Array.isArray(e) ? e : e.bones
	},
	getBoneByName: function(e, t) {
		for (var i = 0, n = this.getBones(t); i < n.length; i++)
			if (e === n[i].name) return n[i]
	},
	getNearestBone: function(e, t) {
		for (; e.isBone;) {
			if (-1 !== t.indexOf(e.name)) return e;
			e = e.parent
		}
	},
	findBoneTrackData: function(e, t) {
		for (var i = /\[(.*)\]\.(.*)/, n = {
				name: e
			}, r = 0; r < t.length; ++r) {
			var a = i.exec(t[r].name);
			a && e === a[1] && (n[a[2]] = r)
		}
		return n
	},
	getEqualsBonesNames: function(e, t) {
		var i = this.getBones(e),
			n = this.getBones(t),
			r = [];
		e: for (var a = 0; a < i.length; a++)
			for (var s = i[a].name, o = 0; o < n.length; o++)
				if (s === n[o].name) {
					r.push(s);
					continue e
				}
		return r
	},
	clone: function(e) {
		var t = new Map,
			i = new Map,
			n = e.clone();
		return parallelTraverse(e, n, (function(e, n) {
			t.set(n, e), i.set(e, n)
		})), n.traverse((function(e) {
			if (e.isSkinnedMesh) {
				var n = e,
					r = t.get(e),
					a = r.skeleton.bones;
				n.skeleton = r.skeleton.clone(), n.bindMatrix.copy(r.bindMatrix), n.skeleton.bones = a.map((function(e) {
					return i.get(e)
				})), n.bind(n.skeleton, n.bindMatrix)
			}
		})), n
	}
};

function parallelTraverse(e, t, i) {
	i(e, t);
	for (var n = 0; n < e.children.length; n++) parallelTraverse(e.children[n], t.children[n], i)
}
const $prepared = Symbol("prepared"),
	$prepare = Symbol("prepare"),
	$preparedGLTF = Symbol("preparedGLTF"),
	$clone = Symbol("clone");
class GLTFInstance {
	constructor(e) {
		this[$preparedGLTF] = e
	}
	static prepare(e) {
		if (null == e.scene) throw new Error("Model does not have a scene");
		if (e[$prepared]) return e;
		const t = this[$prepare](e);
		return t[$prepared] = !0, t
	}
	static[$prepare](e) {
		const {
			scene: t
		} = e, i = [t];
		return Object.assign(Object.assign({}, e), {
			scene: t,
			scenes: i
		})
	}
	get parser() {
		return this[$preparedGLTF].parser
	}
	get animations() {
		return this[$preparedGLTF].animations
	}
	get scene() {
		return this[$preparedGLTF].scene
	}
	get scenes() {
		return this[$preparedGLTF].scenes
	}
	get cameras() {
		return this[$preparedGLTF].cameras
	}
	get asset() {
		return this[$preparedGLTF].asset
	}
	get userData() {
		return this[$preparedGLTF].userData
	}
	clone() {
		return new(0, this.constructor)(this[$clone]())
	}
	dispose() {
		this.scenes.forEach(e => {
			e.traverse(e => {
				if (!e.isMesh) return;
				const t = e;
				(Array.isArray(t.material) ? t.material : [t.material]).forEach(e => {
					e.dispose()
				}), t.geometry.dispose()
			})
		})
	} [$clone]() {
		const e = this[$preparedGLTF],
			t = SkeletonUtils.clone(this.scene),
			i = [t],
			n = e.userData ? Object.assign({}, e.userData) : {};
		return Object.assign(Object.assign({}, e), {
			scene: t,
			scenes: i,
			userData: n
		})
	}
}
const alphaChunk = "\n#ifdef ALPHATEST\n\n    if ( diffuseColor.a < ALPHATEST ) discard;\n    diffuseColor.a = 1.0;\n\n#endif\n",
	$threeGLTF = Symbol("threeGLTF"),
	$gltf = Symbol("gltf"),
	$gltfElementMap = Symbol("gltfElementMap"),
	$threeObjectMap = Symbol("threeObjectMap"),
	$parallelTraverseThreeScene = Symbol("parallelTraverseThreeScene"),
	$correlateOriginalThreeGLTF = Symbol("correlateOriginalThreeGLTF"),
	$correlateCloneThreeGLTF = Symbol("correlateCloneThreeGLTF");
class CorrelatedSceneGraph {
	constructor(e, t, i, n) {
		this[$threeGLTF] = e, this[$gltf] = t, this[$gltfElementMap] = n, this[$threeObjectMap] = i
	}
	static from(e, t) {
		return null != t ? this[$correlateCloneThreeGLTF](e, t) : this[$correlateOriginalThreeGLTF](e)
	}
	static[$correlateOriginalThreeGLTF](e) {
		const t = e.parser.json,
			{
				associations: i
			} = e.parser,
			n = new Map,
			r = {
				name: "Default"
			},
			a = {
				type: "materials",
				index: -1
			};
		return i.forEach((e, i) => {
			null == e && (a.index < 0 && (null == t.materials && (t.materials = []), a.index = t.materials.length, t.materials.push(r)), e = a);
			const {
				type: s,
				index: o
			} = e, l = (t[s] || [])[o];
			if (null == l) return;
			let c = n.get(l);
			null == c && (c = new Set, n.set(l, c)), c.add(i)
		}), new CorrelatedSceneGraph(e, t, i, n)
	}
	static[$correlateCloneThreeGLTF](e, t) {
		const i = t.threeGLTF,
			n = t.gltf,
			r = JSON.parse(JSON.stringify(n)),
			a = new Map,
			s = new Map;
		for (let n = 0; n < i.scenes.length; n++) this[$parallelTraverseThreeScene](i.scenes[n], e.scenes[n], (e, i) => {
			const n = t.threeObjectMap.get(e);
			if (null != n) {
				const {
					type: e,
					index: t
				} = n, o = r[e][t];
				a.set(i, {
					type: e,
					index: t
				});
				const l = s.get(o) || new Set;
				l.add(i), s.set(o, l)
			}
		});
		return new CorrelatedSceneGraph(e, r, a, s)
	}
	static[$parallelTraverseThreeScene](e, t, i) {
		const n = (e, t) => {
			if (i(e, t), e.isObject3D) {
				if (e.isMesh)
					if (Array.isArray(e.material))
						for (let i = 0; i < e.material.length; ++i) n(e.material[i], t.material[i]);
					else n(e.material, t.material);
				for (let i = 0; i < e.children.length; ++i) n(e.children[i], t.children[i])
			}
		};
		n(e, t)
	}
	get threeGLTF() {
		return this[$threeGLTF]
	}
	get gltf() {
		return this[$gltf]
	}
	get gltfElementMap() {
		return this[$gltfElementMap]
	}
	get threeObjectMap() {
		return this[$threeObjectMap]
	}
	loadVariant(e, t = (() => {})) {
		const i = new Set;
		return this.threeGLTF.scene.traverse(async n => {
			const {
				gltfExtensions: r
			} = n.userData;
			if (!n.isMesh || null == r) return;
			const a = r.KHR_materials_variants;
			if (null == a) return;
			let s = -1;
			for (const t of a.mappings)
				if (t.variants.indexOf(e) >= 0) {
					s = t.material;
					break
				} if (s < 0) return;
			const o = await this.threeGLTF.parser.getDependency("material", s);
			i.add(s), n.material = o, this.threeGLTF.parser.assignFinalMaterial(n), t();
			const l = this.gltf.materials[s];
			let c = this.gltfElementMap.get(l);
			null == c && (c = new Set, this.gltfElementMap.set(l, c)), c.add(n.material)
		}), i
	}
}
const $cloneAndPatchMaterial = Symbol("cloneAndPatchMaterial"),
	$correlatedSceneGraph = Symbol("correlatedSceneGraph");
class ModelViewerGLTFInstance extends GLTFInstance {
	static[$prepare](e) {
		const t = super[$prepare](e);
		null == t[$correlatedSceneGraph] && (t[$correlatedSceneGraph] = CorrelatedSceneGraph.from(t));
		const {
			scene: i
		} = t, n = [];
		i.traverse(e => {
			if (e.renderOrder = 1e3, e.frustumCulled = !1, e.name || (e.name = e.uuid), !e.isMesh) return;
			e.castShadow = !0;
			const t = e;
			let i = !1;
			(Array.isArray(t.material) ? t.material : [t.material]).forEach(e => {
				e.isMeshStandardMaterial && (e.transparent && 2 === e.side && (i = !0, e.side = 0), Renderer.singleton.roughnessMipmapper.generateMipmaps(e))
			}), i && n.push(t)
		});
		for (const e of n) {
			const t = (Array.isArray(e.material) ? e.material : [e.material]).map(e => {
					const t = e.clone();
					return t.side = 1, t
				}),
				i = Array.isArray(e.material) ? t : t[0],
				n = new Mesh(e.geometry, i);
			n.renderOrder = -1, e.add(n)
		}
		return t
	}
	get correlatedSceneGraph() {
		return this[$preparedGLTF][$correlatedSceneGraph]
	} [$clone]() {
		const e = super[$clone](),
			t = new Map;
		return e.scene.traverse(e => {
			if (e.isMesh) {
				const i = e;
				Array.isArray(i.material) ? i.material = i.material.map(e => this[$cloneAndPatchMaterial](e, t)) : null != i.material && (i.material = this[$cloneAndPatchMaterial](i.material, t))
			}
		}), e[$correlatedSceneGraph] = CorrelatedSceneGraph.from(e, this.correlatedSceneGraph), e
	} [$cloneAndPatchMaterial](e, t) {
		var i, n;
		if (t.has(e.uuid)) return t.get(e.uuid);
		const r = e.clone();
		null != e.map && (r.map = e.map.clone(), r.map.needsUpdate = !0), null != e.normalMap && (r.normalMap = e.normalMap.clone(), r.normalMap.needsUpdate = !0), null != e.emissiveMap && (r.emissiveMap = e.emissiveMap.clone(), r.emissiveMap.needsUpdate = !0), e.isGLTFSpecularGlossinessMaterial ? (null != e.specularMap && (r.specularMap = null === (i = e.specularMap) || void 0 === i ? void 0 : i.clone(), r.specularMap.needsUpdate = !0), null != e.glossinessMap && (r.glossinessMap = null === (n = e.glossinessMap) || void 0 === n ? void 0 : n.clone(), r.glossinessMap.needsUpdate = !0)) : (e.metalnessMap === e.aoMap ? r.metalnessMap = r.aoMap : null != e.metalnessMap && (r.metalnessMap = e.metalnessMap.clone(), r.metalnessMap.needsUpdate = !0), e.roughnessMap === e.aoMap ? r.roughnessMap = r.aoMap : e.roughnessMap === e.metalnessMap ? r.roughnessMap = r.metalnessMap : null != e.roughnessMap && (r.roughnessMap = e.roughnessMap.clone(), r.roughnessMap.needsUpdate = !0));
		const a = e.onBeforeCompile;
		return r.onBeforeCompile = e.isGLTFSpecularGlossinessMaterial ? e => {
			a(e, void 0), e.fragmentShader = e.fragmentShader.replace("#include <alphatest_fragment>", alphaChunk)
		} : e => {
			e.fragmentShader = e.fragmentShader.replace("#include <alphatest_fragment>", alphaChunk), a(e, void 0)
		}, r.shadowSide = 0, r.transparent && (r.depthWrite = !1), r.alphaTest || r.transparent || (r.alphaTest = -.5), t.set(e.uuid, r), r
	}
}
var RGBELoader = function(e) {
	DataTextureLoader.call(this, e), this.type = 1009
};
RGBELoader.prototype = Object.assign(Object.create(DataTextureLoader.prototype), {
	constructor: RGBELoader,
	parse: function(e) {
		var t = function(e, t) {
				switch (e) {
					case 1:
						console.error("THREE.RGBELoader Read Error: " + (t || ""));
						break;
					case 2:
						console.error("THREE.RGBELoader Write Error: " + (t || ""));
						break;
					case 3:
						console.error("THREE.RGBELoader Bad File Format: " + (t || ""));
						break;
					default:
					case 4:
						console.error("THREE.RGBELoader: Error: " + (t || ""))
				}
				return -1
			},
			i = function(e, t, i) {
				t = t || 1024;
				for (var n = e.pos, r = -1, a = 0, s = "", o = String.fromCharCode.apply(null, new Uint16Array(e.subarray(n, n + 128))); 0 > (r = o.indexOf("\n")) && a < t && n < e.byteLength;) s += o, a += o.length, n += 128, o += String.fromCharCode.apply(null, new Uint16Array(e.subarray(n, n + 128)));
				return -1 < r && (!1 !== i && (e.pos += a + r + 1), s + o.slice(0, r))
			},
			n = function(e, t, i, n) {
				var r = e[t + 3],
					a = Math.pow(2, r - 128) / 255;
				i[n + 0] = e[t + 0] * a, i[n + 1] = e[t + 1] * a, i[n + 2] = e[t + 2] * a
			},
			r = function(e, t, i, n) {
				var r = e[t + 3],
					a = Math.pow(2, r - 128) / 255;
				i[n + 0] = DataUtils.toHalfFloat(e[t + 0] * a), i[n + 1] = DataUtils.toHalfFloat(e[t + 1] * a), i[n + 2] = DataUtils.toHalfFloat(e[t + 2] * a)
			},
			a = new Uint8Array(e);
		a.pos = 0;
		var s = function(e) {
			var n, r, a = /^\s*GAMMA\s*=\s*(\d+(\.\d+)?)\s*$/,
				s = /^\s*EXPOSURE\s*=\s*(\d+(\.\d+)?)\s*$/,
				o = /^\s*FORMAT=(\S+)\s*$/,
				l = /^\s*\-Y\s+(\d+)\s+\+X\s+(\d+)\s*$/,
				c = {
					valid: 0,
					string: "",
					comments: "",
					programtype: "RGBE",
					format: "",
					gamma: 1,
					exposure: 1,
					width: 0,
					height: 0
				};
			if (e.pos >= e.byteLength || !(n = i(e))) return t(1, "no header found");
			if (!(r = n.match(/^#\?(\S+)/))) return t(3, "bad initial token");
			for (c.valid |= 1, c.programtype = r[1], c.string += n + "\n"; !1 !== (n = i(e));)
				if (c.string += n + "\n", "#" !== n.charAt(0)) {
					if ((r = n.match(a)) && (c.gamma = parseFloat(r[1], 10)), (r = n.match(s)) && (c.exposure = parseFloat(r[1], 10)), (r = n.match(o)) && (c.valid |= 2, c.format = r[1]), (r = n.match(l)) && (c.valid |= 4, c.height = parseInt(r[1], 10), c.width = parseInt(r[2], 10)), 2 & c.valid && 4 & c.valid) break
				} else c.comments += n + "\n";
			return 2 & c.valid ? 4 & c.valid ? c : t(3, "missing image size specifier") : t(3, "missing format specifier")
		}(a);
		if (-1 !== s) {
			var o = s.width,
				l = s.height,
				c = function(e, i, n) {
					var r, a, s, o, l, c, h, u, d, p, m, A, g, f = i,
						v = n;
					if (f < 8 || f > 32767 || 2 !== e[0] || 2 !== e[1] || 128 & e[2]) return new Uint8Array(e);
					if (f !== (e[2] << 8 | e[3])) return t(3, "wrong scanline width");
					if (!(r = new Uint8Array(4 * i * n)).length) return t(4, "unable to allocate buffer space");
					for (a = 0, s = 0, u = 4 * f, g = new Uint8Array(4), c = new Uint8Array(u); v > 0 && s < e.byteLength;) {
						if (s + 4 > e.byteLength) return t(1);
						if (g[0] = e[s++], g[1] = e[s++], g[2] = e[s++], g[3] = e[s++], 2 != g[0] || 2 != g[1] || (g[2] << 8 | g[3]) != f) return t(3, "bad rgbe scanline format");
						for (h = 0; h < u && s < e.byteLength;) {
							if ((A = (o = e[s++]) > 128) && (o -= 128), 0 === o || h + o > u) return t(3, "bad scanline data");
							if (A)
								for (l = e[s++], d = 0; d < o; d++) c[h++] = l;
							else c.set(e.subarray(s, s + o), h), h += o, s += o
						}
						for (p = f, d = 0; d < p; d++) m = 0, r[a] = c[d + m], m += f, r[a + 1] = c[d + m], m += f, r[a + 2] = c[d + m], m += f, r[a + 3] = c[d + m], a += 4;
						v--
					}
					return r
				}(a.subarray(a.pos), o, l);
			if (-1 !== c) {
				switch (this.type) {
					case 1009:
						var h = c,
							u = 1023,
							d = 1009;
						break;
					case 1015:
						for (var p = c.length / 4 * 3, m = new Float32Array(p), A = 0; A < p; A++) n(c, 4 * A, m, 3 * A);
						h = m, u = 1022, d = 1015;
						break;
					case 1016:
						p = c.length / 4 * 3;
						var g = new Uint16Array(p);
						for (A = 0; A < p; A++) r(c, 4 * A, g, 3 * A);
						h = g, u = 1022, d = 1016;
						break;
					default:
						console.error("THREE.RGBELoader: unsupported type: ", this.type)
				}
				return {
					width: o,
					height: l,
					data: h,
					header: s.string,
					gamma: s.gamma,
					exposure: s.exposure,
					format: u,
					type: d
				}
			}
		}
		return null
	},
	setDataType: function(e) {
		return this.type = e, this
	},
	load: function(e, t, i, n) {
		return DataTextureLoader.prototype.load.call(this, e, (function(e, i) {
			switch (e.type) {
				case 1009:
					e.encoding = 3002, e.minFilter = 1003, e.magFilter = 1003, e.generateMipmaps = !1, e.flipY = !0;
					break;
				case 1015:
				case 1016:
					e.encoding = 3e3, e.minFilter = 1006, e.magFilter = 1006, e.generateMipmaps = !1, e.flipY = !0
			}
			t && t(e, i)
		}), i, n)
	}
});
class EnvironmentScene extends Scene {
	constructor() {
		super(), this.position.y = -3.5;
		const e = new BoxGeometry;
		e.deleteAttribute("uv");
		const t = new MeshStandardMaterial({
				metalness: 0,
				side: 1
			}),
			i = new MeshStandardMaterial({
				metalness: 0
			}),
			n = new PointLight(16777215, 500, 28, 2);
		n.position.set(.418, 16.199, .3), this.add(n);
		const r = new Mesh(e, t);
		r.position.set(-.757, 13.219, .717), r.scale.set(31.713, 28.305, 28.591), this.add(r);
		const a = new Mesh(e, i);
		a.position.set(-10.906, 2.009, 1.846), a.rotation.set(0, -.195, 0), a.scale.set(2.328, 7.905, 4.651), this.add(a);
		const s = new Mesh(e, i);
		s.position.set(-5.607, -.754, -.758), s.rotation.set(0, .994, 0), s.scale.set(1.97, 1.534, 3.955), this.add(s);
		const o = new Mesh(e, i);
		o.position.set(6.167, .857, 7.803), o.rotation.set(0, .561, 0), o.scale.set(3.927, 6.285, 3.687), this.add(o);
		const l = new Mesh(e, i);
		l.position.set(-2.017, .018, 6.124), l.rotation.set(0, .333, 0), l.scale.set(2.002, 4.566, 2.064), this.add(l);
		const c = new Mesh(e, i);
		c.position.set(2.291, -.756, -2.621), c.rotation.set(0, -.286, 0), c.scale.set(1.546, 1.552, 1.496), this.add(c);
		const h = new Mesh(e, i);
		h.position.set(-2.193, -.369, -5.547), h.rotation.set(0, .516, 0), h.scale.set(3.875, 3.487, 2.986), this.add(h);
		const u = new Mesh(e, this.createAreaLightMaterial(50));
		u.position.set(-16.116, 14.37, 8.208), u.scale.set(.1, 2.428, 2.739), this.add(u);
		const d = new Mesh(e, this.createAreaLightMaterial(50));
		d.position.set(-16.109, 18.021, -8.207), d.scale.set(.1, 2.425, 2.751), this.add(d);
		const p = new Mesh(e, this.createAreaLightMaterial(17));
		p.position.set(14.904, 12.198, -1.832), p.scale.set(.15, 4.265, 6.331), this.add(p);
		const m = new Mesh(e, this.createAreaLightMaterial(43));
		m.position.set(-.462, 8.89, 14.52), m.scale.set(4.38, 5.441, .088), this.add(m);
		const A = new Mesh(e, this.createAreaLightMaterial(20));
		A.position.set(3.235, 11.486, -12.541), A.scale.set(2.5, 2, .1), this.add(A);
		const g = new Mesh(e, this.createAreaLightMaterial(100));
		g.position.set(0, 20, 0), g.scale.set(1, .1, 1), this.add(g)
	}
	createAreaLightMaterial(e) {
		const t = new MeshBasicMaterial;
		return t.color.setScalar(e), t
	}
}
class EnvironmentSceneAlt extends Scene {
	constructor() {
		super(), this.position.y = -3.5;
		const e = new BoxGeometry;
		e.deleteAttribute("uv");
		const t = new MeshStandardMaterial({
				metalness: 0,
				side: 1
			}),
			i = new MeshStandardMaterial({
				metalness: 0
			}),
			n = new PointLight(16777215, 500, 28, 2);
		n.position.set(.5, 16.5, .5), this.add(n);
		const r = new Mesh(e, t);
		r.position.set(0, 13.2, 0), r.scale.set(31.5, 28.5, 31.5), this.add(r);
		const a = new Mesh(e, i);
		a.position.set(-10.906, -1, 1.846), a.rotation.set(0, -.195, 0), a.scale.set(2.328, 7.905, 4.651), this.add(a);
		const s = new Mesh(e, i);
		s.position.set(-5.607, -.754, -.758), s.rotation.set(0, .994, 0), s.scale.set(1.97, 1.534, 3.955), this.add(s);
		const o = new Mesh(e, i);
		o.position.set(6.167, -.16, 7.803), o.rotation.set(0, .561, 0), o.scale.set(3.927, 6.285, 3.687), this.add(o);
		const l = new Mesh(e, i);
		l.position.set(-2.017, .018, 6.124), l.rotation.set(0, .333, 0), l.scale.set(2.002, 4.566, 2.064), this.add(l);
		const c = new Mesh(e, i);
		c.position.set(2.291, -.756, -2.621), c.rotation.set(0, -.286, 0), c.scale.set(1.546, 1.552, 1.496), this.add(c);
		const h = new Mesh(e, i);
		h.position.set(-2.193, -.369, -5.547), h.rotation.set(0, .516, 0), h.scale.set(3.875, 3.487, 2.986), this.add(h);
		const u = new Mesh(e, this.createAreaLightMaterial(18));
		u.position.set(-14, 9, .1), u.scale.set(.1, 5, 5), this.add(u);
		const d = new Mesh(e, this.createAreaLightMaterial(18));
		d.position.set(14, 9, .1), d.scale.set(.1, 5, 5), this.add(d);
		const p = new Mesh(e, this.createAreaLightMaterial(18));
		p.position.set(0, 9, 14), p.scale.set(5, 5, .1), this.add(p);
		const m = new Mesh(e, this.createAreaLightMaterial(18));
		m.position.set(0, 9, -14), m.scale.set(5, 5, .1), this.add(m);
		const A = new Mesh(e, this.createAreaLightMaterial(18));
		A.position.set(0, 20, 0), A.scale.set(1.5, .1, 1.5), this.add(A)
	}
	createAreaLightMaterial(e) {
		const t = new MeshBasicMaterial;
		return t.color.setScalar(e), t
	}
}
const GENERATED_SIGMA = .04,
	HDR_FILE_RE = /\.hdr(\.js)?$/,
	ldrLoader = new TextureLoader,
	hdrLoader = new RGBELoader,
	userData = {
		url: null
	};
class TextureUtils extends EventDispatcher {
	constructor(e) {
		super(), this.generatedEnvironmentMap = null, this.generatedEnvironmentMapAlt = null, this.skyboxCache = new Map, this.environmentMapCache = new Map, this.PMREMGenerator = new PMREMGenerator(e)
	}
	async load(e, t = (() => {})) {
		try {
			const i = HDR_FILE_RE.test(e),
				n = i ? hdrLoader : ldrLoader,
				r = await new Promise((i, r) => n.load(e, i, e => {
					t(e.loaded / e.total * .9)
				}, r));
			return t(1), this.addMetadata(r, e), r.mapping = 303, i ? (r.encoding = 3002, r.minFilter = 1003, r.magFilter = 1003, r.flipY = !0) : r.encoding = 3007, r
		} finally {
			t && t(1)
		}
	}
	async generateEnvironmentMapAndSkybox(e = null, t = null, i = {}) {
		const {
			progressTracker: n
		} = i, r = null != n ? n.beginActivity() : () => {}, a = "neutral" === t;
		!0 === a && (t = null);
		const s = deserializeUrl(t);
		try {
			let t, i = Promise.resolve(null);
			e && (i = this.loadSkyboxFromUrl(e, n)), t = s ? this.loadEnvironmentMapFromUrl(s, n) : e ? this.loadEnvironmentMapFromUrl(e, n) : !0 === a ? this.loadGeneratedEnvironmentMapAlt() : this.loadGeneratedEnvironmentMap();
			let [o, l] = await Promise.all([t, i]);
			if (null == o) throw new Error("Failed to load environment map.");
			return {
				environmentMap: o,
				skybox: l
			}
		} finally {
			r(1)
		}
	}
	addMetadata(e, t) {
		null != e && (e.userData = Object.assign(Object.assign({}, userData), {
			url: t
		}))
	}
	loadSkyboxFromUrl(e, t) {
		if (!this.skyboxCache.has(e)) {
			const i = t ? t.beginActivity() : () => {},
				n = this.load(e, i);
			this.skyboxCache.set(e, n)
		}
		return this.skyboxCache.get(e)
	}
	loadEnvironmentMapFromUrl(e, t) {
		if (!this.environmentMapCache.has(e)) {
			const i = this.loadSkyboxFromUrl(e, t).then(t => {
				const i = this.PMREMGenerator.fromEquirectangular(t);
				return this.addMetadata(i.texture, e), i
			});
			this.PMREMGenerator.compileEquirectangularShader(), this.environmentMapCache.set(e, i)
		}
		return this.environmentMapCache.get(e)
	}
	loadGeneratedEnvironmentMap() {
		if (null == this.generatedEnvironmentMap) {
			const e = new EnvironmentScene;
			this.generatedEnvironmentMap = this.PMREMGenerator.fromScene(e, .04), this.addMetadata(this.generatedEnvironmentMap.texture, null)
		}
		return Promise.resolve(this.generatedEnvironmentMap)
	}
	loadGeneratedEnvironmentMapAlt() {
		if (null == this.generatedEnvironmentMapAlt) {
			const e = new EnvironmentSceneAlt;
			this.generatedEnvironmentMapAlt = this.PMREMGenerator.fromScene(e, .04), this.addMetadata(this.generatedEnvironmentMapAlt.texture, null)
		}
		return Promise.resolve(this.generatedEnvironmentMapAlt)
	}
	async dispose() {
		const e = [];
		this.environmentMapCache.forEach(t => {
			e.push(t)
		}), this.environmentMapCache.clear();
		for (const t of e) try {
			(await t).dispose()
		} catch (e) {}
		null != this.generatedEnvironmentMap && (this.generatedEnvironmentMap.dispose(), this.generatedEnvironmentMap = null), null != this.generatedEnvironmentMapAlt && (this.generatedEnvironmentMapAlt.dispose(), this.generatedEnvironmentMapAlt = null)
	}
}
const DURATION_DECAY = .2,
	LOW_FRAME_DURATION_MS = 18,
	HIGH_FRAME_DURATION_MS = 26,
	MAX_AVG_CHANGE_MS = 2,
	SCALE_STEPS = [1, .79, .62, .5, .4, .31, .25],
	DEFAULT_LAST_STEP = 3;
class Renderer extends EventDispatcher {
	constructor(e) {
		super(), this.loader = new CachingGLTFLoader(ModelViewerGLTFInstance), this.width = 0, this.height = 0, this.dpr = 1, this.debugger = null, this.scenes = new Set, this.multipleScenesVisible = !1, this.scaleStep = 0, this.lastStep = 3, this.avgFrameDuration = 22, this.onWebGLContextLost = e => {
			this.dispatchEvent({
				type: "contextlost",
				sourceEvent: e
			})
		}, this.dpr = resolveDpr(), this.canvasElement = document.createElement("canvas"), this.canvasElement.id = "webgl-canvas", this.canvas3D = this.canvasElement, this.canvas3D.addEventListener("webglcontextlost", this.onWebGLContextLost);
		try {
			this.threeRenderer = new WebGL1Renderer({
				canvas: this.canvas3D,
				alpha: !0,
				antialias: !0,
				powerPreference: "high-performance",
				preserveDrawingBuffer: !0
			}), this.threeRenderer.autoClear = !0, this.threeRenderer.outputEncoding = 3007, this.threeRenderer.physicallyCorrectLights = !0, this.threeRenderer.setPixelRatio(1), this.threeRenderer.shadowMap.enabled = !0, this.threeRenderer.shadowMap.type = 2, this.threeRenderer.shadowMap.autoUpdate = !1, this.debugger = null != e && e.debug ? new Debugger(this) : null, this.threeRenderer.debug = {
				checkShaderErrors: !!this.debugger
			}, this.threeRenderer.toneMapping = 4
		} catch (e) {
			console.warn(e)
		}
		this.arRenderer = new ARRenderer(this), this.textureUtils = this.canRender ? new TextureUtils(this.threeRenderer) : null, this.roughnessMipmapper = new RoughnessMipmapper(this.threeRenderer), CachingGLTFLoader.initializeKTX2Loader(this.threeRenderer), this.updateRendererSize(), this.lastTick = performance.now(), this.avgFrameDuration = 0
	}
	static get singleton() {
		return this._singleton
	}
	static resetSingleton() {
		this._singleton.dispose(), this._singleton = new Renderer({
			debug: isDebugMode()
		})
	}
	get canRender() {
		return null != this.threeRenderer
	}
	get scaleFactor() {
		return SCALE_STEPS[this.scaleStep]
	}
	set minScale(e) {
		let t = 1;
		for (; t < SCALE_STEPS.length && !(SCALE_STEPS[t] < e);) ++t;
		this.lastStep = t - 1
	}
	updateRendererSize() {
		const e = resolveDpr();
		if (e !== this.dpr)
			for (const e of this.scenes) {
				const {
					element: t
				} = e;
				t[$updateSize](t.getBoundingClientRect())
			}
		let t = 0,
			i = 0;
		for (const e of this.scenes) t = Math.max(t, e.width), i = Math.max(i, e.height);
		if (t === this.width && i === this.height && e === this.dpr) return;
		this.width = t, this.height = i, this.dpr = e, this.canRender && this.threeRenderer.setSize(t * e, i * e, !1);
		const n = this.scaleFactor,
			r = t / n,
			a = i / n;
		this.canvasElement.style.width = r + "px", this.canvasElement.style.height = a + "px";
		for (const n of this.scenes) {
			const {
				canvas: s
			} = n;
			s.width = Math.round(t * e), s.height = Math.round(i * e), s.style.width = r + "px", s.style.height = a + "px", n.isDirty = !0
		}
	}
	updateRendererScale() {
		const e = this.scaleStep;
		if (this.avgFrameDuration > 26 && this.scaleStep < this.lastStep ? ++this.scaleStep : this.avgFrameDuration < 18 && this.scaleStep > 0 && --this.scaleStep, e == this.scaleStep) return;
		const t = this.scaleFactor;
		this.avgFrameDuration = 22;
		const i = this.width / t,
			n = this.height / t;
		this.canvasElement.style.width = i + "px", this.canvasElement.style.height = n + "px";
		for (const e of this.scenes) {
			const {
				style: t
			} = e.canvas;
			t.width = i + "px", t.height = n + "px", e.isDirty = !0
		}
	}
	registerScene(e) {
		this.scenes.add(e);
		const {
			canvas: t
		} = e, i = this.scaleFactor;
		t.width = Math.round(this.width * this.dpr), t.height = Math.round(this.height * this.dpr), t.style.width = this.width / i + "px", t.style.height = this.height / i + "px", this.multipleScenesVisible && t.classList.add("show"), e.isDirty = !0, this.canRender && this.scenes.size > 0 && this.threeRenderer.setAnimationLoop(e => this.render(e)), null != this.debugger && this.debugger.addScene(e)
	}
	unregisterScene(e) {
		this.scenes.delete(e), this.canRender && 0 === this.scenes.size && this.threeRenderer.setAnimationLoop(null), null != this.debugger && this.debugger.removeScene(e)
	}
	displayCanvas(e) {
		return this.multipleScenesVisible ? e.element[$canvas] : this.canvasElement
	}
	selectCanvas() {
		let e = 0,
			t = null;
		for (const i of this.scenes) {
			const {
				element: n
			} = i;
			n.modelIsVisible && (++e, t = n[$userInputElement])
		}
		const i = e > 1 || !1,
			{
				canvasElement: n
			} = this;
		if (i !== this.multipleScenesVisible || !i && n.parentElement !== t) {
			this.multipleScenesVisible = i, i && n.classList.remove("show");
			for (const e of this.scenes) {
				const r = e.element[$userInputElement],
					a = e.element[$canvas];
				i ? (a.classList.add("show"), e.isDirty = !0) : r === t && (r.appendChild(n), n.classList.add("show"), a.classList.remove("show"), e.isDirty = !0)
			}
		}
	}
	orderedScenes() {
		const e = [];
		for (const t of [!1, !0])
			for (const i of this.scenes) i.element.modelIsVisible === t && e.push(i);
		return e
	}
	get isPresenting() {
		return this.arRenderer.isPresenting
	}
	preRender(e, t, i) {
		const {
			element: n,
			exposure: r
		} = e;
		n[$tick](t, i);
		const a = "number" == typeof r && !self.isNaN(r);
		this.threeRenderer.toneMappingExposure = a ? r : 1, e.isShadowDirty() && (this.threeRenderer.shadowMap.needsUpdate = !0)
	}
	render(e) {
		const t = e - this.lastTick;
		if (this.lastTick = e, !this.canRender || this.isPresenting) return;
		this.avgFrameDuration += clamp(.2 * (t - this.avgFrameDuration), -2, 2), this.selectCanvas(), this.updateRendererSize(), this.updateRendererScale();
		const {
			dpr: i,
			scaleFactor: n
		} = this;
		for (const r of this.orderedScenes()) {
			if (!r.element[$sceneIsReady]()) continue;
			if (this.preRender(r, e, t), !r.isDirty) continue;
			if (r.isDirty = !1, ++r.renderCount, !r.element.modelIsVisible && !this.multipleScenesVisible)
				for (const e of this.scenes) e.element.modelIsVisible && (e.isDirty = !0);
			const a = Math.min(Math.ceil(r.width * n * i), this.canvas3D.width),
				s = Math.min(Math.ceil(r.height * n * i), this.canvas3D.height);
			if (this.threeRenderer.setRenderTarget(null), this.threeRenderer.setViewport(0, Math.floor(this.height * i) - s, a, s), this.threeRenderer.render(r, r.getCamera()), this.multipleScenesVisible) {
				null == r.context && r.createContext(); {
					const e = r.context;
					e.clearRect(0, 0, a, s), e.drawImage(this.canvas3D, 0, 0, a, s, 0, 0, a, s)
				}
			}
		}
	}
	dispose() {
		null != this.textureUtils && this.textureUtils.dispose(), null != this.threeRenderer && this.threeRenderer.dispose(), this.textureUtils = null, this.threeRenderer = null, this.scenes.clear(), this.canvas3D.removeEventListener("webglcontextlost", this.onWebGLContextLost)
	}
}
Renderer._singleton = new Renderer({
	debug: isDebugMode()
});
const dataUrlToBlob = async e => new Promise((t, i) => {
	const n = e.match(/data:(.*);/);
	if (!n) return i(new Error(e + " is not a valid data Url"));
	const r = n[1],
		a = e.replace(/data:image\/\w+;base64,/, ""),
		s = atob(a),
		o = [];
	for (let e = 0; e < s.length; e += 512) {
		const t = s.slice(e, e + 512),
			i = new Array(t.length);
		for (let e = 0; e < t.length; e++) i[e] = t.charCodeAt(e);
		const n = new Uint8Array(i);
		o.push(n)
	}
	t(new Blob(o, {
		type: r
	}))
});
var _a$8, _b$7;
const $ongoingActivities = Symbol("ongoingActivities"),
	$announceTotalProgress = Symbol("announceTotalProgress"),
	$eventDelegate = Symbol("eventDelegate"),
	ACTIVITY_PROGRESS_WEIGHT = .5;
class ProgressTracker {
	constructor() {
		this[_a$8] = document.createDocumentFragment(), this.addEventListener = (...e) => this[$eventDelegate].addEventListener(...e), this.removeEventListener = (...e) => this[$eventDelegate].removeEventListener(...e), this.dispatchEvent = (...e) => this[$eventDelegate].dispatchEvent(...e), this[_b$7] = new Set
	}
	get ongoingActivityCount() {
		return this[$ongoingActivities].size
	}
	beginActivity() {
		const e = {
			progress: 0
		};
		return this[$ongoingActivities].add(e), 1 === this.ongoingActivityCount && this[$announceTotalProgress](), t => {
			let i;
			return i = Math.max(clamp(t, 0, 1), e.progress), i !== e.progress && (e.progress = i, this[$announceTotalProgress]()), e.progress
		}
	} [(_a$8 = $eventDelegate, _b$7 = $ongoingActivities, $announceTotalProgress)]() {
		let e = 0,
			t = 0,
			i = 0;
		for (const n of this[$ongoingActivities]) {
			const {
				progress: r
			} = n;
			e += r * (.5 / Math.pow(2, t++)), 1 === r && i++
		}
		i === this.ongoingActivityCount && (e = 1, this[$ongoingActivities].clear()), this.dispatchEvent(new CustomEvent("progress", {
			detail: {
				totalProgress: e
			}
		}))
	}
}
var _a$7, _b$6, _c$3, _d$1, _e, _f, _g, _h, _j, _k, __decorate$7 = function(e, t, i, n) {
	for (var r, a = arguments.length, s = a < 3 ? t : null === n ? n = Object.getOwnPropertyDescriptor(t, i) : n, o = e.length - 1; o >= 0; o--)(r = e[o]) && (s = (a < 3 ? r(s) : a > 3 ? r(t, i, s) : r(t, i)) || s);
	return a > 3 && s && Object.defineProperty(t, i, s), s
};
const CLEAR_MODEL_TIMEOUT_MS = 1e3,
	FALLBACK_SIZE_UPDATE_THRESHOLD_MS = 50,
	ANNOUNCE_MODEL_VISIBILITY_DEBOUNCE_THRESHOLD = 0,
	UNSIZED_MEDIA_WIDTH = 300,
	UNSIZED_MEDIA_HEIGHT = 150,
	blobCanvas = document.createElement("canvas");
let blobContext = null;
const $template = Symbol("template"),
	$fallbackResizeHandler = Symbol("fallbackResizeHandler"),
	$defaultAriaLabel = Symbol("defaultAriaLabel"),
	$resizeObserver = Symbol("resizeObserver"),
	$intersectionObserver = Symbol("intersectionObserver"),
	$clearModelTimeout = Symbol("clearModelTimeout"),
	$onContextLost = Symbol("onContextLost"),
	$loaded = Symbol("loaded"),
	$updateSize = Symbol("updateSize"),
	$isElementInViewport = Symbol("isElementInViewport"),
	$announceModelVisibility = Symbol("announceModelVisibility"),
	$ariaLabel = Symbol("ariaLabel"),
	$loadedTime = Symbol("loadedTime"),
	$updateSource = Symbol("updateSource"),
	$markLoaded = Symbol("markLoaded"),
	$container = Symbol("container"),
	$userInputElement = Symbol("input"),
	$canvas = Symbol("canvas"),
	$scene = Symbol("scene"),
	$needsRender = Symbol("needsRender"),
	$tick = Symbol("tick"),
	$onModelLoad = Symbol("onModelLoad"),
	$onResize = Symbol("onResize"),
	$renderer = Symbol("renderer"),
	$progressTracker = Symbol("progressTracker"),
	$getLoaded = Symbol("getLoaded"),
	$getModelIsVisible = Symbol("getModelIsVisible"),
	$shouldAttemptPreload = Symbol("shouldAttemptPreload"),
	$sceneIsReady = Symbol("sceneIsReady"),
	$hasTransitioned = Symbol("hasTransitioned"),
	toVector3D = e => ({
		x: e.x,
		y: e.y,
		z: e.z,
		toString() {
			return `${this.x}m ${this.y}m ${this.z}m`
		}
	});
class ModelViewerElementBase extends UpdatingElement {
	constructor() {
		super(), this.alt = null, this.src = null, this[_a$7] = !1, this[_b$6] = !1, this[_c$3] = 0, this[_d$1] = null, this[_e] = debounce(() => {
			const e = this.getBoundingClientRect();
			this[$updateSize](e)
		}, 50), this[_f] = debounce(e => {
			const t = this.modelIsVisible;
			t !== e && this.dispatchEvent(new CustomEvent("model-visibility", {
				detail: {
					visible: t
				}
			}))
		}, 0), this[_g] = null, this[_h] = null, this[_j] = new ProgressTracker, this[_k] = e => {
			this.dispatchEvent(new CustomEvent("error", {
				detail: {
					type: "webglcontextlost",
					sourceError: e.sourceEvent
				}
			}))
		};
		const e = this.constructor.template;
		window.ShadyCSS && window.ShadyCSS.styleElement(this, {}), this.attachShadow({
			mode: "open"
		});
		const t = this.shadowRoot;
		let i, n;
		if (t.appendChild(e.content.cloneNode(!0)), this[$container] = t.querySelector(".container"), this[$userInputElement] = t.querySelector(".userInput"), this[$canvas] = t.querySelector("canvas"), this[$defaultAriaLabel] = this[$userInputElement].getAttribute("aria-label"), this.isConnected) {
			const e = this.getBoundingClientRect();
			i = e.width, n = e.height
		} else i = 300, n = 150;
		this[$scene] = new ModelScene({
			canvas: this[$canvas],
			element: this,
			width: i,
			height: n
		}), this[$scene].addEventListener("model-load", async e => {
			this[$markLoaded](), this[$onModelLoad](), await timePasses(), this.dispatchEvent(new CustomEvent("load", {
				detail: {
					url: e.url
				}
			}))
		}), Promise.resolve().then(() => {
			this[$updateSize](this.getBoundingClientRect())
		}), HAS_RESIZE_OBSERVER && (this[$resizeObserver] = new ResizeObserver(e => {
			if (!this[$renderer].isPresenting)
				for (let t of e) t.target === this && this[$updateSize](t.contentRect)
		})), HAS_INTERSECTION_OBSERVER ? this[$intersectionObserver] = new IntersectionObserver(e => {
			for (let t of e)
				if (t.target === this) {
					const e = this.modelIsVisible;
					this[$isElementInViewport] = t.isIntersecting, this[$announceModelVisibility](e), this[$isElementInViewport] && !this[$sceneIsReady]() && this[$updateSource]()
				}
		}, {
			root: null,
			rootMargin: "0px",
			threshold: 0
		}) : this[$isElementInViewport] = !0
	}
	static get is() {
		return "model-viewer"
	}
	static get template() {
		return this.hasOwnProperty($template) || (this[$template] = makeTemplate(this.is)), this[$template]
	}
	static set modelCacheSize(e) {
		CachingGLTFLoader[$evictionPolicy].evictionThreshold = e
	}
	static get modelCacheSize() {
		return CachingGLTFLoader[$evictionPolicy].evictionThreshold
	}
	static set minimumRenderScale(e) {
		e > 1 && console.warn("<model-viewer> minimumRenderScale has been clamped to a maximum value of 1."), e <= 0 && console.warn("<model-viewer> minimumRenderScale has been clamped to a minimum value of 0.25."), Renderer.singleton.minScale = e
	}
	static get minimumRenderScale() {
		return Renderer.singleton.minScale
	}
	get loaded() {
		return this[$getLoaded]()
	}
	get[(_a$7 = $isElementInViewport, _b$6 = $loaded, _c$3 = $loadedTime, _d$1 = $clearModelTimeout, _e = $fallbackResizeHandler, _f = $announceModelVisibility, _g = $resizeObserver, _h = $intersectionObserver, _j = $progressTracker, $renderer)]() {
		return Renderer.singleton
	}
	get modelIsVisible() {
		return this[$getModelIsVisible]()
	}
	connectedCallback() {
		super.connectedCallback && super.connectedCallback(), HAS_RESIZE_OBSERVER ? this[$resizeObserver].observe(this) : self.addEventListener("resize", this[$fallbackResizeHandler]), HAS_INTERSECTION_OBSERVER && this[$intersectionObserver].observe(this);
		const e = this[$renderer];
		e.addEventListener("contextlost", this[$onContextLost]), e.registerScene(this[$scene]), null != this[$clearModelTimeout] && (self.clearTimeout(this[$clearModelTimeout]), this[$clearModelTimeout] = null, this.requestUpdate("src", null))
	}
	disconnectedCallback() {
		super.disconnectedCallback && super.disconnectedCallback(), HAS_RESIZE_OBSERVER ? this[$resizeObserver].unobserve(this) : self.removeEventListener("resize", this[$fallbackResizeHandler]), HAS_INTERSECTION_OBSERVER && this[$intersectionObserver].unobserve(this);
		const e = this[$renderer];
		e.removeEventListener("contextlost", this[$onContextLost]), e.unregisterScene(this[$scene]), this[$clearModelTimeout] = self.setTimeout(() => {
			this[$scene].reset()
		}, 1e3)
	}
	updated(e) {
		if (super.updated(e), e.has("src") && (null == this.src ? (this[$loaded] = !1, this[$loadedTime] = 0, this[$scene].reset()) : this.src !== this[$scene].url && (this[$loaded] = !1, this[$loadedTime] = 0, this[$updateSource]())), e.has("alt")) {
			const e = null == this.alt ? this[$defaultAriaLabel] : this.alt;
			this[$userInputElement].setAttribute("aria-label", e)
		}
	}
	toDataURL(e, t) {
		return this[$renderer].displayCanvas(this[$scene]).toDataURL(e, t)
	}
	async toBlob(e) {
		const t = e ? e.mimeType : void 0,
			i = e ? e.qualityArgument : void 0,
			n = e ? e.idealAspect : void 0,
			{
				width: r,
				height: a,
				fieldOfViewAspect: s,
				aspect: o
			} = this[$scene],
			{
				dpr: l,
				scaleFactor: c
			} = this[$renderer];
		let h = r * c * l,
			u = a * c * l,
			d = 0,
			p = 0;
		if (!0 === n)
			if (s > o) {
				const e = u;
				u = Math.round(h / s), p = (e - u) / 2
			} else {
				const e = h;
				h = Math.round(u * s), d = (e - h) / 2
			} blobCanvas.width = h, blobCanvas.height = u;
		try {
			return new Promise(async (e, n) => (null == blobContext && (blobContext = blobCanvas.getContext("2d")), blobContext.drawImage(this[$renderer].displayCanvas(this[$scene]), d, p, h, u, 0, 0, h, u), !blobCanvas.msToBlob || t && "image/png" !== t ? blobCanvas.toBlob ? void blobCanvas.toBlob(t => {
				if (!t) return n(new Error("Unable to retrieve canvas blob"));
				e(t)
			}, t, i) : e(await dataUrlToBlob(blobCanvas.toDataURL(t, i))) : e(blobCanvas.msToBlob())))
		} finally {
			this[$updateSize]({
				width: r,
				height: a
			})
		}
	}
	get[$ariaLabel]() {
		return null == this.alt || "null" === this.alt ? this[$defaultAriaLabel] : this.alt
	} [$getLoaded]() {
		return this[$loaded]
	} [$getModelIsVisible]() {
		return this.loaded && this[$isElementInViewport]
	} [$hasTransitioned]() {
		return this.modelIsVisible
	} [$shouldAttemptPreload]() {
		return !!this.src && this[$isElementInViewport]
	} [$sceneIsReady]() {
		return this[$loaded]
	} [$updateSize]({
		width: e,
		height: t
	}) {
		this[$container].style.width = e + "px", this[$container].style.height = t + "px", this[$onResize]({
			width: parseFloat(e),
			height: parseFloat(t)
		})
	} [$tick](e, t) {} [$markLoaded]() {
		this[$loaded] || (this[$loaded] = !0, this[$loadedTime] = performance.now())
	} [$needsRender]() {
		this[$scene].isDirty = !0
	} [$onModelLoad]() {} [$onResize](e) {
		this[$scene].setSize(e.width, e.height)
	}
	async [(_k = $onContextLost, $updateSource)]() {
		if (this.loaded || !this[$shouldAttemptPreload]()) return;
		const e = this[$progressTracker].beginActivity(),
			t = this.src;
		try {
			await this[$scene].setSource(t, t => e(.8 * t));
			const i = {
				url: t
			};
			this.dispatchEvent(new CustomEvent("preload", {
				detail: i
			}))
		} catch (e) {
			this.dispatchEvent(new CustomEvent("error", {
				detail: e
			}))
		} finally {
			e(.9), requestAnimationFrame(() => {
				requestAnimationFrame(() => {
					e(1)
				})
			})
		}
	}
}
__decorate$7([property({
	type: String
})], ModelViewerElementBase.prototype, "alt", void 0), __decorate$7([property({
	type: String
})], ModelViewerElementBase.prototype, "src", void 0);
var __decorate$6 = function(e, t, i, n) {
	for (var r, a = arguments.length, s = a < 3 ? t : null === n ? n = Object.getOwnPropertyDescriptor(t, i) : n, o = e.length - 1; o >= 0; o--)(r = e[o]) && (s = (a < 3 ? r(s) : a > 3 ? r(t, i, s) : r(t, i)) || s);
	return a > 3 && s && Object.defineProperty(t, i, s), s
};
const MILLISECONDS_PER_SECOND = 1e3,
	$changeAnimation = Symbol("changeAnimation"),
	$paused = Symbol("paused"),
	AnimationMixin = e => {
		var t;
		class i extends e {
			constructor() {
				super(...arguments), this.autoplay = !1, this.animationName = void 0, this.animationCrossfadeDuration = 300, this[t] = !0
			}
			get availableAnimations() {
				return this.loaded ? this[$scene].animationNames : []
			}
			get duration() {
				return this[$scene].duration
			}
			get paused() {
				return this[$paused]
			}
			get currentTime() {
				return this[$scene].animationTime
			}
			set currentTime(e) {
				this[$scene].animationTime = e, this[$renderer].threeRenderer.shadowMap.needsUpdate = !0, this[$needsRender]()
			}
			pause() {
				this[$paused] || (this[$paused] = !0, this[$renderer].threeRenderer.shadowMap.autoUpdate = !1, this.dispatchEvent(new CustomEvent("pause")))
			}
			play() {
				this[$paused] && this.availableAnimations.length > 0 && (this[$paused] = !1, this[$renderer].threeRenderer.shadowMap.autoUpdate = !0, this[$scene].hasActiveAnimation || this[$changeAnimation](), this.dispatchEvent(new CustomEvent("play")))
			} [(t = $paused, $onModelLoad)]() {
				super[$onModelLoad](), this[$paused] = !0, this.autoplay && (this[$changeAnimation](), this.play())
			} [$tick](e, t) {
				super[$tick](e, t), !this[$paused] && this[$hasTransitioned]() && (this[$scene].updateAnimation(t / 1e3), this[$needsRender]())
			}
			updated(e) {
				super.updated(e), e.has("autoplay") && this.autoplay && this.play(), e.has("animationName") && this[$changeAnimation]()
			}
			async [$updateSource]() {
				return this[$scene].stopAnimation(), super[$updateSource]()
			} [$changeAnimation]() {
				this[$scene].playAnimation(this.animationName, this.animationCrossfadeDuration / 1e3), this[$paused] && (this[$scene].updateAnimation(0), this[$needsRender]())
			}
		}
		return __decorate$6([property({
			type: Boolean
		})], i.prototype, "autoplay", void 0), __decorate$6([property({
			type: String,
			attribute: "animation-name"
		})], i.prototype, "animationName", void 0), __decorate$6([property({
			type: Number,
			attribute: "animation-crossfade-duration"
		})], i.prototype, "animationCrossfadeDuration", void 0), i
	},
	$annotationRenderer = Symbol("annotationRenderer"),
	$hotspotMap = Symbol("hotspotMap"),
	$mutationCallback = Symbol("mutationCallback"),
	$observer = Symbol("observer"),
	$addHotspot = Symbol("addHotspot"),
	$removeHotspot = Symbol("removeHotspot"),
	pixelPosition = new Vector2,
	worldToModel = new Matrix4,
	worldToModelNormal = new Matrix3,
	AnnotationMixin = e => {
		var t, i, n, r;
		class a extends e {
			constructor(...e) {
				super(...e), this[t] = new CSS2DRenderer, this[i] = new Map, this[n] = e => {
					e.forEach(e => {
						e instanceof MutationRecord && "childList" !== e.type || (e.addedNodes.forEach(e => {
							this[$addHotspot](e)
						}), e.removedNodes.forEach(e => {
							this[$removeHotspot](e)
						}), this[$needsRender]())
					})
				}, this[r] = new MutationObserver(this[$mutationCallback]);
				const {
					domElement: a
				} = this[$annotationRenderer], {
					style: s
				} = a;
				s.display = "none", s.pointerEvents = "none", s.position = "absolute", s.top = "0", this.shadowRoot.querySelector(".default").appendChild(a)
			}
			connectedCallback() {
				super.connectedCallback();
				for (let e = 0; e < this.children.length; ++e) this[$addHotspot](this.children[e]);
				const {
					ShadyDOM: e
				} = self;
				null == e ? this[$observer].observe(this, {
					childList: !0
				}) : this[$observer] = e.observeChildren(this, this[$mutationCallback])
			}
			disconnectedCallback() {
				super.disconnectedCallback();
				const {
					ShadyDOM: e
				} = self;
				null == e ? this[$observer].disconnect() : e.unobserveChildren(this[$observer])
			}
			updateHotspot(e) {
				const t = this[$hotspotMap].get(e.name);
				null != t && (t.updatePosition(e.position), t.updateNormal(e.normal), this[$needsRender]())
			}
			positionAndNormalFromPoint(e, t) {
				const i = this[$scene],
					{
						width: n,
						height: r,
						target: a
					} = i;
				pixelPosition.set(e / n, t / r).multiplyScalar(2).subScalar(1), pixelPosition.y *= -1;
				const s = i.positionAndNormalFromPoint(pixelPosition);
				if (null == s) return null;
				worldToModel.copy(a.matrixWorld).invert();
				const o = toVector3D(s.position.applyMatrix4(worldToModel));
				worldToModelNormal.getNormalMatrix(worldToModel);
				return {
					position: o,
					normal: toVector3D(s.normal.applyNormalMatrix(worldToModelNormal))
				}
			} [(t = $annotationRenderer, i = $hotspotMap, n = $mutationCallback, r = $observer, $tick)](e, t) {
				super[$tick](e, t);
				const i = this[$scene],
					n = i.getCamera();
				i.isDirty && (i.updateHotspots(n.position), this[$annotationRenderer].domElement.style.display = "", this[$annotationRenderer].render(i, n))
			} [$onResize](e) {
				super[$onResize](e), this[$annotationRenderer].setSize(e.width, e.height)
			} [$addHotspot](e) {
				if (!(e instanceof HTMLElement && 0 === e.slot.indexOf("hotspot"))) return;
				let t = this[$hotspotMap].get(e.slot);
				null != t ? t.increment() : (t = new Hotspot({
					name: e.slot,
					position: e.dataset.position,
					normal: e.dataset.normal
				}), this[$hotspotMap].set(e.slot, t), this[$scene].addHotspot(t), this[$annotationRenderer].domElement.appendChild(t.element)), this[$scene].isDirty = !0
			} [$removeHotspot](e) {
				if (!(e instanceof HTMLElement)) return;
				const t = this[$hotspotMap].get(e.slot);
				t && (t.decrement() && (this[$scene].removeHotspot(t), this[$hotspotMap].delete(e.slot)), this[$scene].isDirty = !0)
			}
		}
		return a
	},
	enumerationDeserializer = e => t => {
		try {
			const i = parseExpressions(t),
				n = (i.length ? i[0].terms : []).filter(e => e && "ident" === e.type).map(e => e.value).filter(t => e.indexOf(t) > -1),
				r = new Set;
			for (const e of n) r.add(e);
			return r
		} catch (e) {}
		return new Set
	};
var __decorate$5 = function(e, t, i, n) {
	for (var r, a = arguments.length, s = a < 3 ? t : null === n ? n = Object.getOwnPropertyDescriptor(t, i) : n, o = e.length - 1; o >= 0; o--)(r = e[o]) && (s = (a < 3 ? r(s) : a > 3 ? r(t, i, s) : r(t, i)) || s);
	return a > 3 && s && Object.defineProperty(t, i, s), s
};
let isWebXRBlocked = !1,
	isSceneViewerBlocked = !1;
const noArViewerSigil = "#model-viewer-no-ar-fallback",
	deserializeARModes = enumerationDeserializer(["quick-look", "scene-viewer", "webxr", "none"]),
	DEFAULT_AR_MODES = "webxr scene-viewer quick-look",
	ARMode = {
		QUICK_LOOK: "quick-look",
		SCENE_VIEWER: "scene-viewer",
		WEBXR: "webxr",
		NONE: "none"
	},
	$arButtonContainer = Symbol("arButtonContainer"),
	$enterARWithWebXR = Symbol("enterARWithWebXR"),
	$openSceneViewer = Symbol("openSceneViewer"),
	$openIOSARQuickLook = Symbol("openIOSARQuickLook"),
	$canActivateAR = Symbol("canActivateAR"),
	$arMode = Symbol("arMode"),
	$arModes = Symbol("arModes"),
	$arAnchor = Symbol("arAnchor"),
	$preload = Symbol("preload"),
	$onARButtonContainerClick = Symbol("onARButtonContainerClick"),
	$onARStatus = Symbol("onARStatus"),
	$onARTap = Symbol("onARTap"),
	$selectARMode = Symbol("selectARMode"),
	ARMixin = e => {
		var t, i, n, r, a, s, o, l, c;
		class h extends e {
			constructor() {
				super(...arguments), this.ar = !1, this.arScale = "auto", this.arPlacement = "floor", this.arModes = DEFAULT_AR_MODES, this.iosSrc = null, this[t] = !1, this[i] = this.shadowRoot.querySelector(".ar-button"), this[n] = document.createElement("a"), this[r] = new Set, this[a] = ARMode.NONE, this[s] = !1, this[o] = e => {
					e.preventDefault(), this.activateAR()
				}, this[l] = ({
					status: e
				}) => {
					e !== ARStatus.NOT_PRESENTING && this[$renderer].arRenderer.presentedScene !== this[$scene] || (this.setAttribute("ar-status", e), this.dispatchEvent(new CustomEvent("ar-status", {
						detail: {
							status: e
						}
					})))
				}, this[c] = e => {
					"_apple_ar_quicklook_button_tapped" == e.data && this.dispatchEvent(new CustomEvent("quick-look-button-tapped"))
				}
			}
			get canActivateAR() {
				return this[$arMode] !== ARMode.NONE
			}
			connectedCallback() {
				super.connectedCallback(), this[$renderer].arRenderer.addEventListener("status", this[$onARStatus]), this.setAttribute("ar-status", ARStatus.NOT_PRESENTING), this[$arAnchor].addEventListener("message", this[$onARTap])
			}
			disconnectedCallback() {
				super.disconnectedCallback(), this[$renderer].arRenderer.removeEventListener("status", this[$onARStatus]), this[$arAnchor].removeEventListener("message", this[$onARTap])
			}
			async update(e) {
				super.update(e), e.has("arScale") && (this[$scene].canScale = "fixed" !== this.arScale), e.has("arPlacement") && (this[$scene].setShadowIntensity(this[$scene].shadowIntensity), this[$needsRender]()), (e.has("ar") || e.has("arModes") || e.has("iosSrc")) && (e.has("arModes") && (this[$arModes] = deserializeARModes(this.arModes)), this[$selectARMode]())
			}
			async activateAR() {
				switch (this[$arMode]) {
					case ARMode.QUICK_LOOK:
						this[$openIOSARQuickLook]();
						break;
					case ARMode.WEBXR:
						await this[$enterARWithWebXR]();
						break;
					case ARMode.SCENE_VIEWER:
						this[$openSceneViewer]();
						break;
					default:
						console.warn("No AR Mode can be activated. This is probably due to missing configuration or device capabilities")
				}
			}
			async [(t = $canActivateAR, i = $arButtonContainer, n = $arAnchor, r = $arModes, a = $arMode, s = $preload, o = $onARButtonContainerClick, l = $onARStatus, c = $onARTap, $selectARMode)]() {
				if (this[$arMode] = ARMode.NONE, this.ar) {
					const e = [];
					this[$arModes].forEach(t => {
						e.push(t)
					});
					for (const t of e) {
						if ("webxr" === t && IS_WEBXR_AR_CANDIDATE && !isWebXRBlocked && await this[$renderer].arRenderer.supportsPresentation()) {
							this[$arMode] = ARMode.WEBXR;
							break
						}
						if ("scene-viewer" === t && IS_SCENEVIEWER_CANDIDATE && !isSceneViewerBlocked) {
							this[$arMode] = ARMode.SCENE_VIEWER;
							break
						}
						if ("quick-look" === t && this.iosSrc && IS_AR_QUICKLOOK_CANDIDATE) {
							this[$arMode] = ARMode.QUICK_LOOK;
							break
						}
					}
				}
				if (this.canActivateAR) this[$arButtonContainer].classList.add("enabled"), this[$arButtonContainer].addEventListener("click", this[$onARButtonContainerClick]);
				else if (this[$arButtonContainer].classList.contains("enabled")) {
					this[$arButtonContainer].removeEventListener("click", this[$onARButtonContainerClick]), this[$arButtonContainer].classList.remove("enabled");
					const e = ARStatus.FAILED;
					this.setAttribute("ar-status", e), this.dispatchEvent(new CustomEvent("ar-status", {
						detail: {
							status: e
						}
					}))
				}
			}
			async [$enterARWithWebXR]() {
				console.log("Attempting to present in AR..."), this[$loaded] || (this[$preload] = !0, this[$updateSource](), await waitForEvent(this, "load"), this[$preload] = !1);
				try {
					this[$arButtonContainer].removeEventListener("click", this[$onARButtonContainerClick]);
					const {
						arRenderer: e
					} = this[$renderer];
					e.placeOnWall = "wall" === this.arPlacement, await e.present(this[$scene])
				} catch (e) {
					console.warn("Error while trying to present to AR"), console.error(e), await this[$renderer].arRenderer.stopPresenting(), isWebXRBlocked = !0, await this[$selectARMode](), this.activateAR()
				} finally {
					this[$selectARMode]()
				}
			} [$shouldAttemptPreload]() {
				return super[$shouldAttemptPreload]() || this[$preload]
			} [$openSceneViewer]() {
				const e = self.location.toString(),
					t = new URL(e),
					i = new URL(this.src, e),
					n = new URLSearchParams(i.search);
				if (t.hash = noArViewerSigil, n.set("mode", "ar_only"), n.has("disable_occlusion") || n.set("disable_occlusion", "true"), "fixed" === this.arScale && n.set("resizable", "false"), "wall" === this.arPlacement && n.set("enable_vertical_placement", "true"), n.has("sound")) {
					const t = new URL(n.get("sound"), e);
					n.set("sound", t.toString())
				}
				if (n.has("link")) {
					const t = new URL(n.get("link"), e);
					n.set("link", t.toString())
				}
				const r = `intent://arvr.google.com/scene-viewer/1.0?${n.toString()+"&file="+encodeURIComponent(i.toString())}#Intent;scheme=https;package=com.google.ar.core;action=android.intent.action.VIEW;S.browser_fallback_url=${encodeURIComponent(t.toString())};end;`;
				self.addEventListener("hashchange", () => {
					self.location.hash === noArViewerSigil && (isSceneViewerBlocked = !0, self.history.back(), this[$selectARMode]())
				}, {
					once: !0
				}), this[$arAnchor].setAttribute("href", r), this[$arAnchor].click()
			} [$openIOSARQuickLook]() {
				const e = new URL(this.iosSrc, self.location.toString());
				"fixed" === this.arScale && (e.hash && (e.hash += "&"), e.hash += "allowsContentScaling=0");
				const t = this[$arAnchor];
				t.setAttribute("rel", "ar");
				const i = document.createElement("img");
				t.appendChild(i), t.setAttribute("href", e.toString()), t.click(), t.removeChild(i)
			}
		}
		return __decorate$5([property({
			type: Boolean,
			attribute: "ar"
		})], h.prototype, "ar", void 0), __decorate$5([property({
			type: String,
			attribute: "ar-scale"
		})], h.prototype, "arScale", void 0), __decorate$5([property({
			type: String,
			attribute: "ar-placement"
		})], h.prototype, "arPlacement", void 0), __decorate$5([property({
			type: String,
			attribute: "ar-modes"
		})], h.prototype, "arModes", void 0), __decorate$5([property({
			type: String,
			attribute: "ios-src"
		})], h.prototype, "iosSrc", void 0), h
	};
var _a$6, _b$5, _c$2;
const $evaluate = Symbol("evaluate"),
	$lastValue = Symbol("lastValue");
class Evaluator {
	constructor() {
		this[_a$6] = null
	}
	static evaluatableFor(e, t = ZERO) {
		if (e instanceof Evaluator) return e;
		if ("number" === e.type) return "%" === e.unit ? new PercentageEvaluator(e, t) : e;
		switch (e.name.value) {
			case "calc":
				return new CalcEvaluator(e, t);
			case "env":
				return new EnvEvaluator(e)
		}
		return ZERO
	}
	static evaluate(e) {
		return e instanceof Evaluator ? e.evaluate() : e
	}
	static isConstant(e) {
		return !(e instanceof Evaluator) || e.isConstant
	}
	static applyIntrinsics(e, t) {
		const {
			basis: i,
			keywords: n
		} = t, {
			auto: r
		} = n;
		return i.map((t, i) => {
			const a = null == r[i] ? t : r[i];
			let s = e[i] ? e[i] : a;
			if ("ident" === s.type) {
				const e = s.value;
				e in n && (s = n[e][i])
			}
			return null != s && "ident" !== s.type || (s = a), "%" === s.unit ? numberNode(s.number / 100 * t.number, t.unit) : (s = normalizeUnit(s, t), s.unit !== t.unit ? t : s)
		})
	}
	get isConstant() {
		return !1
	}
	evaluate() {
		return this.isConstant && null != this[$lastValue] || (this[$lastValue] = this[$evaluate]()), this[$lastValue]
	}
}
_a$6 = $lastValue;
const $percentage = Symbol("percentage"),
	$basis = Symbol("basis");
class PercentageEvaluator extends Evaluator {
	constructor(e, t) {
		super(), this[$percentage] = e, this[$basis] = t
	}
	get isConstant() {
		return !0
	} [$evaluate]() {
		return numberNode(this[$percentage].number / 100 * this[$basis].number, this[$basis].unit)
	}
}
const $identNode = Symbol("identNode");
class EnvEvaluator extends Evaluator {
	constructor(e) {
		super(), this[_b$5] = null;
		const t = e.arguments.length ? e.arguments[0].terms[0] : null;
		null != t && "ident" === t.type && (this[$identNode] = t)
	}
	get isConstant() {
		return !1
	} [(_b$5 = $identNode, $evaluate)]() {
		if (null != this[$identNode]) switch (this[$identNode].value) {
			case "window-scroll-y":
				return {
					type: "number", number: window.pageYOffset / (Math.max(document.body.scrollHeight, document.body.offsetHeight, document.documentElement.clientHeight, document.documentElement.scrollHeight, document.documentElement.offsetHeight) - window.innerHeight) || 0, unit: null
				}
		}
		return ZERO
	}
}
const IS_MULTIPLICATION_RE = /[\*\/]/,
	$evaluator = Symbol("evalutor");
class CalcEvaluator extends Evaluator {
	constructor(e, t = ZERO) {
		if (super(), this[_c$2] = null, 1 !== e.arguments.length) return;
		const i = e.arguments[0].terms.slice(),
			n = [];
		for (; i.length;) {
			const e = i.shift();
			if (n.length > 0) {
				const i = n[n.length - 1];
				if ("operator" === i.type && IS_MULTIPLICATION_RE.test(i.value)) {
					const i = n.pop(),
						r = n.pop();
					if (null == r) return;
					n.push(new OperatorEvaluator(i, Evaluator.evaluatableFor(r, t), Evaluator.evaluatableFor(e, t)));
					continue
				}
			}
			n.push("operator" === e.type ? e : Evaluator.evaluatableFor(e, t))
		}
		for (; n.length > 2;) {
			const [e, i, r] = n.splice(0, 3);
			if ("operator" !== i.type) return;
			n.unshift(new OperatorEvaluator(i, Evaluator.evaluatableFor(e, t), Evaluator.evaluatableFor(r, t)))
		}
		1 === n.length && (this[$evaluator] = n[0])
	}
	get isConstant() {
		return null == this[$evaluator] || Evaluator.isConstant(this[$evaluator])
	} [(_c$2 = $evaluator, $evaluate)]() {
		return null != this[$evaluator] ? Evaluator.evaluate(this[$evaluator]) : ZERO
	}
}
const $operator = Symbol("operator"),
	$left = Symbol("left"),
	$right = Symbol("right");
class OperatorEvaluator extends Evaluator {
	constructor(e, t, i) {
		super(), this[$operator] = e, this[$left] = t, this[$right] = i
	}
	get isConstant() {
		return Evaluator.isConstant(this[$left]) && Evaluator.isConstant(this[$right])
	} [$evaluate]() {
		const e = normalizeUnit(Evaluator.evaluate(this[$left])),
			t = normalizeUnit(Evaluator.evaluate(this[$right])),
			{
				number: i,
				unit: n
			} = e,
			{
				number: r,
				unit: a
			} = t;
		if (null != a && null != n && a != n) return ZERO;
		const s = n || a;
		let o;
		switch (this[$operator].value) {
			case "+":
				o = i + r;
				break;
			case "-":
				o = i - r;
				break;
			case "/":
				o = i / r;
				break;
			case "*":
				o = i * r;
				break;
			default:
				return ZERO
		}
		return {
			type: "number",
			number: o,
			unit: s
		}
	}
}
const $evaluatables = Symbol("evaluatables"),
	$intrinsics = Symbol("intrinsics");
class StyleEvaluator extends Evaluator {
	constructor(e, t) {
		super(), this[$intrinsics] = t;
		const i = e[0],
			n = null != i ? i.terms : [];
		this[$evaluatables] = t.basis.map((e, t) => {
			const i = n[t];
			return null == i ? {
				type: "ident",
				value: "auto"
			} : "ident" === i.type ? i : Evaluator.evaluatableFor(i, e)
		})
	}
	get isConstant() {
		for (const e of this[$evaluatables])
			if (!Evaluator.isConstant(e)) return !1;
		return !0
	} [$evaluate]() {
		const e = this[$evaluatables].map(e => Evaluator.evaluate(e));
		return Evaluator.applyIntrinsics(e, this[$intrinsics]).map(e => e.number)
	}
}
var _a$5, _b$4, _c$1, _d;
const $instances = Symbol("instances"),
	$activateListener = Symbol("activateListener"),
	$deactivateListener = Symbol("deactivateListener"),
	$notifyInstances = Symbol("notifyInstances"),
	$notify = Symbol("notify"),
	$scrollCallback = Symbol("callback");
class ScrollObserver {
	constructor(e) {
		this[$scrollCallback] = e
	}
	static[$notifyInstances]() {
		for (const e of ScrollObserver[$instances]) e[$notify]()
	}
	static[(_a$5 = $instances, $activateListener)]() {
		window.addEventListener("scroll", this[$notifyInstances], {
			passive: !0
		})
	}
	static[$deactivateListener]() {
		window.removeEventListener("scroll", this[$notifyInstances])
	}
	observe() {
		0 === ScrollObserver[$instances].size && ScrollObserver[$activateListener](), ScrollObserver[$instances].add(this)
	}
	disconnect() {
		ScrollObserver[$instances].delete(this), 0 === ScrollObserver[$instances].size && ScrollObserver[$deactivateListener]()
	} [$notify]() {
		this[$scrollCallback]()
	}
}
ScrollObserver[_a$5] = new Set;
const $computeStyleCallback = Symbol("computeStyleCallback"),
	$astWalker = Symbol("astWalker"),
	$dependencies = Symbol("dependencies"),
	$onScroll = Symbol("onScroll");
class StyleEffector {
	constructor(e) {
		this[_b$4] = {}, this[_c$1] = new ASTWalker(["function"]), this[_d] = () => {
			this[$computeStyleCallback]({
				relatedState: "window-scroll"
			})
		}, this[$computeStyleCallback] = e
	}
	observeEffectsFor(e) {
		const t = {},
			i = this[$dependencies];
		this[$astWalker].walk(e, e => {
			const {
				name: n
			} = e, r = e.arguments[0].terms[0];
			if ("env" === n.value && null != r && "ident" === r.type) switch (r.value) {
				case "window-scroll-y":
					if (null == t["window-scroll"]) {
						const e = "window-scroll" in i ? i["window-scroll"] : new ScrollObserver(this[$onScroll]);
						e.observe(), delete i["window-scroll"], t["window-scroll"] = e
					}
			}
		});
		for (const e in i) {
			i[e].disconnect()
		}
		this[$dependencies] = t
	}
	dispose() {
		for (const e in this[$dependencies]) {
			this[$dependencies][e].disconnect()
		}
	}
}
_b$4 = $dependencies, _c$1 = $astWalker, _d = $onScroll;
const style = e => {
		const t = e.observeEffects || !1,
			i = e.intrinsics instanceof Function ? e.intrinsics : () => e.intrinsics;
		return (n, r) => {
			const a = n.updated,
				s = n.connectedCallback,
				o = n.disconnectedCallback,
				l = Symbol(r + "StyleEffector"),
				c = Symbol(r + "StyleEvaluator"),
				h = Symbol(r + "UpdateEvaluator"),
				u = Symbol(r + "EvaluateAndSync");
			Object.defineProperties(n, {
				[l]: {
					value: null,
					writable: !0
				},
				[c]: {
					value: null,
					writable: !0
				},
				[h]: {
					value: function() {
						const e = parseExpressions(this[r]);
						this[c] = new StyleEvaluator(e, i(this)), null == this[l] && t && (this[l] = new StyleEffector(() => this[u]())), null != this[l] && this[l].observeEffectsFor(e)
					}
				},
				[u]: {
					value: function() {
						if (null == this[c]) return;
						const t = this[c].evaluate();
						this[e.updateHandler](t)
					}
				},
				updated: {
					value: function(e) {
						e.has(r) && (this[h](), this[u]()), a.call(this, e)
					}
				},
				connectedCallback: {
					value: function() {
						s.call(this), this.requestUpdate(r, this[r])
					}
				},
				disconnectedCallback: {
					value: function() {
						o.call(this), null != this[l] && (this[l].dispose(), this[l] = null)
					}
				}
			})
		}
	},
	DEFAULT_OPTIONS = Object.freeze({
		minimumRadius: 0,
		maximumRadius: 1 / 0,
		minimumPolarAngle: Math.PI / 8,
		maximumPolarAngle: Math.PI - Math.PI / 8,
		minimumAzimuthalAngle: -1 / 0,
		maximumAzimuthalAngle: 1 / 0,
		minimumFieldOfView: 10,
		maximumFieldOfView: 45,
		interactionPolicy: "always-allow",
		touchAction: "pan-y"
	}),
	TOUCH_EVENT_RE = /^touch(start|end|move)$/,
	KEYBOARD_ORBIT_INCREMENT = Math.PI / 8,
	ZOOM_SENSITIVITY = .04,
	KeyCode = {
		PAGE_UP: 33,
		PAGE_DOWN: 34,
		LEFT: 37,
		UP: 38,
		RIGHT: 39,
		DOWN: 40
	},
	ChangeSource = {
		USER_INTERACTION: "user-interaction",
		NONE: "none"
	};
class SmoothControls extends EventDispatcher {
	constructor(e, t) {
		super(), this.camera = e, this.element = t, this.sensitivity = 1, this._interactionEnabled = !1, this._disableZoom = !1, this.isUserChange = !1, this.isUserPointing = !1, this.spherical = new Spherical, this.goalSpherical = new Spherical, this.thetaDamper = new Damper, this.phiDamper = new Damper, this.radiusDamper = new Damper, this.logFov = Math.log(DEFAULT_OPTIONS.maximumFieldOfView), this.goalLogFov = this.logFov, this.fovDamper = new Damper, this.pointerIsDown = !1, this.lastPointerPosition = {
			clientX: 0,
			clientY: 0
		}, this.touchMode = "rotate", this.touchDecided = !1, this.onPointerMove = e => {
			if (this.pointerIsDown && this.canInteract) {
				if (TOUCH_EVENT_RE.test(e.type)) {
					const {
						touches: t
					} = e;
					switch (this.touchMode) {
						case "zoom":
							if (this.lastTouches.length > 1 && t.length > 1) {
								const e = .04 * (this.twoTouchDistance(this.lastTouches[0], this.lastTouches[1]) - this.twoTouchDistance(t[0], t[1])) / 10;
								this.userAdjustOrbit(0, 0, e)
							}
							break;
						case "rotate":
							const {
								touchAction: e
							} = this._options;
							if (!this.touchDecided && "none" !== e) {
								this.touchDecided = !0;
								const {
									clientX: i,
									clientY: n
								} = t[0], r = Math.abs(i - this.lastPointerPosition.clientX), a = Math.abs(n - this.lastPointerPosition.clientY);
								if ("pan-y" === e && a > r || "pan-x" === e && r > a) return void(this.touchMode = "scroll")
							}
							this.handleSinglePointerMove(t[0]);
							break;
						case "scroll":
							return
					}
					this.lastTouches = t
				} else this.handleSinglePointerMove(e);
				e.cancelable && e.preventDefault()
			}
		}, this.onPointerDown = e => {
			if (this.pointerIsDown = !0, this.isUserPointing = !1, TOUCH_EVENT_RE.test(e.type)) {
				const {
					touches: t
				} = e;
				switch (this.touchDecided = !1, t.length) {
					default:
					case 1:
						this.touchMode = "rotate", this.handleSinglePointerDown(t[0]);
						break;
					case 2:
						this.touchMode = this._disableZoom ? "scroll" : "zoom"
				}
				this.lastTouches = t
			} else this.handleSinglePointerDown(e)
		}, this.onPointerUp = e => {
			this.element.style.cursor = "grab", this.pointerIsDown = !1, this.isUserPointing && this.dispatchEvent({
				type: "pointer-change-end",
				pointer: Object.assign({}, this.lastPointerPosition)
			})
		}, this.onWheel = e => {
			if (!this.canInteract) return;
			const t = e.deltaY * (1 == e.deltaMode ? 18 : 1) * .04 / 30;
			this.userAdjustOrbit(0, 0, t), e.cancelable && e.preventDefault()
		}, this.onKeyDown = e => {
			let t = !1;
			switch (e.keyCode) {
				case KeyCode.PAGE_UP:
					t = !0, this.userAdjustOrbit(0, 0, .04);
					break;
				case KeyCode.PAGE_DOWN:
					t = !0, this.userAdjustOrbit(0, 0, -.04);
					break;
				case KeyCode.UP:
					t = !0, this.userAdjustOrbit(0, -KEYBOARD_ORBIT_INCREMENT, 0);
					break;
				case KeyCode.DOWN:
					t = !0, this.userAdjustOrbit(0, KEYBOARD_ORBIT_INCREMENT, 0);
					break;
				case KeyCode.LEFT:
					t = !0, this.userAdjustOrbit(-KEYBOARD_ORBIT_INCREMENT, 0, 0);
					break;
				case KeyCode.RIGHT:
					t = !0, this.userAdjustOrbit(KEYBOARD_ORBIT_INCREMENT, 0, 0)
			}
			t && e.cancelable && e.preventDefault()
		}, this._options = Object.assign({}, DEFAULT_OPTIONS), this.setOrbit(0, Math.PI / 2, 1), this.setFieldOfView(100), this.jumpToGoal()
	}
	get interactionEnabled() {
		return this._interactionEnabled
	}
	enableInteraction() {
		if (!1 === this._interactionEnabled) {
			const {
				element: e
			} = this;
			e.addEventListener("mousemove", this.onPointerMove), e.addEventListener("mousedown", this.onPointerDown), this._disableZoom || e.addEventListener("wheel", this.onWheel), e.addEventListener("keydown", this.onKeyDown), e.addEventListener("touchstart", this.onPointerDown, {
				passive: !0
			}), e.addEventListener("touchmove", this.onPointerMove), self.addEventListener("mouseup", this.onPointerUp), self.addEventListener("touchend", this.onPointerUp), this.element.style.cursor = "grab", this._interactionEnabled = !0
		}
	}
	disableInteraction() {
		if (!0 === this._interactionEnabled) {
			const {
				element: e
			} = this;
			e.removeEventListener("mousemove", this.onPointerMove), e.removeEventListener("mousedown", this.onPointerDown), this._disableZoom || e.removeEventListener("wheel", this.onWheel), e.removeEventListener("keydown", this.onKeyDown), e.removeEventListener("touchstart", this.onPointerDown), e.removeEventListener("touchmove", this.onPointerMove), self.removeEventListener("mouseup", this.onPointerUp), self.removeEventListener("touchend", this.onPointerUp), e.style.cursor = "", this._interactionEnabled = !1
		}
	}
	get options() {
		return this._options
	}
	set disableZoom(e) {
		this._disableZoom != e && (this._disableZoom = e, !0 === e ? this.element.removeEventListener("wheel", this.onWheel) : this.element.addEventListener("wheel", this.onWheel))
	}
	getCameraSpherical(e = new Spherical) {
		return e.copy(this.spherical)
	}
	getFieldOfView() {
		return this.camera.fov
	}
	applyOptions(e) {
		Object.assign(this._options, e), this.setOrbit(), this.setFieldOfView(Math.exp(this.goalLogFov))
	}
	updateNearFar(e, t) {
		this.camera.near = Math.max(e, t / 1e3), this.camera.far = t, this.camera.updateProjectionMatrix()
	}
	updateAspect(e) {
		this.camera.aspect = e, this.camera.updateProjectionMatrix()
	}
	setOrbit(e = this.goalSpherical.theta, t = this.goalSpherical.phi, i = this.goalSpherical.radius) {
		const {
			minimumAzimuthalAngle: n,
			maximumAzimuthalAngle: r,
			minimumPolarAngle: a,
			maximumPolarAngle: s,
			minimumRadius: o,
			maximumRadius: l
		} = this._options, {
			theta: c,
			phi: h,
			radius: u
		} = this.goalSpherical, d = clamp(e, n, r);
		isFinite(n) || isFinite(r) || (this.spherical.theta = this.wrapAngle(this.spherical.theta - d) + d);
		const p = clamp(t, a, s),
			m = clamp(i, o, l);
		return (d !== c || p !== h || m !== u) && (this.goalSpherical.theta = d, this.goalSpherical.phi = p, this.goalSpherical.radius = m, this.goalSpherical.makeSafe(), this.isUserChange = !1, !0)
	}
	setRadius(e) {
		this.goalSpherical.radius = e, this.setOrbit()
	}
	setFieldOfView(e) {
		const {
			minimumFieldOfView: t,
			maximumFieldOfView: i
		} = this._options;
		e = clamp(e, t, i), this.goalLogFov = Math.log(e)
	}
	adjustOrbit(e, t, i) {
		const {
			theta: n,
			phi: r,
			radius: a
		} = this.goalSpherical, {
			minimumRadius: s,
			maximumRadius: o,
			minimumFieldOfView: l,
			maximumFieldOfView: c
		} = this._options, h = this.spherical.theta - n, u = Math.PI - .001, d = n - clamp(e, -u - h, u - h), p = r - t, m = 0 === i ? 0 : i > 0 ? (o - a) / (Math.log(c) - this.goalLogFov) : (a - s) / (this.goalLogFov - Math.log(l)), A = a + i * Math.min(isFinite(m) ? m : 1 / 0, o - s);
		if (this.setOrbit(d, p, A), 0 !== i) {
			const e = this.goalLogFov + i;
			this.setFieldOfView(Math.exp(e))
		}
	}
	jumpToGoal() {
		this.update(0, 1e4)
	}
	update(e, t) {
		if (this.isStationary()) return;
		const {
			maximumPolarAngle: i,
			maximumRadius: n
		} = this._options, r = this.spherical.theta - this.goalSpherical.theta;
		Math.abs(r) > Math.PI && !isFinite(this._options.minimumAzimuthalAngle) && !isFinite(this._options.maximumAzimuthalAngle) && (this.spherical.theta -= 2 * Math.sign(r) * Math.PI), this.spherical.theta = this.thetaDamper.update(this.spherical.theta, this.goalSpherical.theta, t, Math.PI), this.spherical.phi = this.phiDamper.update(this.spherical.phi, this.goalSpherical.phi, t, i), this.spherical.radius = this.radiusDamper.update(this.spherical.radius, this.goalSpherical.radius, t, n), this.logFov = this.fovDamper.update(this.logFov, this.goalLogFov, t, 1), this.moveCamera()
	}
	isStationary() {
		return this.goalSpherical.theta === this.spherical.theta && this.goalSpherical.phi === this.spherical.phi && this.goalSpherical.radius === this.spherical.radius && this.goalLogFov === this.logFov
	}
	moveCamera() {
		this.spherical.makeSafe(), this.camera.position.setFromSpherical(this.spherical), this.camera.setRotationFromEuler(new Euler(this.spherical.phi - Math.PI / 2, this.spherical.theta, 0, "YXZ")), this.camera.fov !== Math.exp(this.logFov) && (this.camera.fov = Math.exp(this.logFov), this.camera.updateProjectionMatrix());
		const e = this.isUserChange ? ChangeSource.USER_INTERACTION : ChangeSource.NONE;
		this.dispatchEvent({
			type: "change",
			source: e
		})
	}
	get canInteract() {
		if ("allow-when-focused" == this._options.interactionPolicy) {
			return this.element.getRootNode().activeElement === this.element
		}
		return "always-allow" === this._options.interactionPolicy
	}
	userAdjustOrbit(e, t, i) {
		this.adjustOrbit(e * this.sensitivity, t * this.sensitivity, i), this.isUserChange = !0, this.dispatchEvent({
			type: "change",
			source: ChangeSource.USER_INTERACTION
		})
	}
	wrapAngle(e) {
		const t = (e + Math.PI) / (2 * Math.PI);
		return 2 * (t - Math.floor(t)) * Math.PI - Math.PI
	}
	pixelLengthToSphericalAngle(e) {
		return 2 * Math.PI * e / this.element.clientHeight
	}
	twoTouchDistance(e, t) {
		const {
			clientX: i,
			clientY: n
		} = e, {
			clientX: r,
			clientY: a
		} = t, s = r - i, o = a - n;
		return Math.sqrt(s * s + o * o)
	}
	handleSinglePointerMove(e) {
		const {
			clientX: t,
			clientY: i
		} = e, n = this.pixelLengthToSphericalAngle(t - this.lastPointerPosition.clientX), r = this.pixelLengthToSphericalAngle(i - this.lastPointerPosition.clientY);
		this.lastPointerPosition.clientX = t, this.lastPointerPosition.clientY = i, !1 === this.isUserPointing && (this.isUserPointing = !0, this.dispatchEvent({
			type: "pointer-change-start",
			pointer: Object.assign({}, e)
		})), this.userAdjustOrbit(n, r, 0)
	}
	handleSinglePointerDown(e) {
		this.lastPointerPosition.clientX = e.clientX, this.lastPointerPosition.clientY = e.clientY, this.element.style.cursor = "grabbing"
	}
}
const easeInOutQuad = e => e < .5 ? 2 * e * e : (4 - 2 * e) * e - 1,
	interpolate = (e, t, i = easeInOutQuad) => n => e + (t - e) * i(n),
	sequence = (e, t) => {
		const i = t.reduce((e, t) => e + t, 0),
			n = t.map(e => e / i);
		return t => {
			let i = 0,
				r = 1 / 0,
				a = () => 0;
			for (let s = 0; s < n.length && (r = n[s], a = e[s], !(t <= i + r)); ++s) i += r;
			return a((t - i) / r)
		}
	},
	timeline = (e, t) => {
		const i = [],
			n = [];
		let r = e;
		for (let e = 0; e < t.length; ++e) {
			const a = t[e],
				{
					value: s,
					frames: o
				} = a,
				l = a.ease || easeInOutQuad,
				c = interpolate(r, s, l);
			i.push(c), n.push(o), r = s
		}
		return sequence(i, n)
	};
var __decorate$4 = function(e, t, i, n) {
	for (var r, a = arguments.length, s = a < 3 ? t : null === n ? n = Object.getOwnPropertyDescriptor(t, i) : n, o = e.length - 1; o >= 0; o--)(r = e[o]) && (s = (a < 3 ? r(s) : a > 3 ? r(t, i, s) : r(t, i)) || s);
	return a > 3 && s && Object.defineProperty(t, i, s), s
};
const PROMPT_ANIMATION_TIME = 5e3,
	wiggle = timeline(0, [{
		frames: 5,
		value: -1
	}, {
		frames: 1,
		value: -1
	}, {
		frames: 8,
		value: 1
	}, {
		frames: 1,
		value: 1
	}, {
		frames: 5,
		value: 0
	}, {
		frames: 18,
		value: 0
	}]),
	fade = timeline(0, [{
		frames: 1,
		value: 1
	}, {
		frames: 5,
		value: 1
	}, {
		frames: 1,
		value: 0
	}, {
		frames: 6,
		value: 0
	}]),
	DEFAULT_CAMERA_ORBIT = "0deg 75deg 105%",
	DEFAULT_CAMERA_TARGET = "auto auto auto",
	DEFAULT_FIELD_OF_VIEW = "auto",
	MINIMUM_RADIUS_RATIO = 1.1 * SAFE_RADIUS_RATIO,
	AZIMUTHAL_QUADRANT_LABELS = ["front", "right", "back", "left"],
	POLAR_TRIENT_LABELS = ["upper-", "", "lower-"],
	DEFAULT_INTERACTION_PROMPT_THRESHOLD = 3e3,
	INTERACTION_PROMPT = "Use mouse, touch or arrow keys to control the camera!",
	InteractionPromptStrategy = {
		AUTO: "auto",
		WHEN_FOCUSED: "when-focused",
		NONE: "none"
	},
	InteractionPromptStyle = {
		BASIC: "basic",
		WIGGLE: "wiggle"
	},
	InteractionPolicy = {
		ALWAYS_ALLOW: "always-allow",
		WHEN_FOCUSED: "allow-when-focused"
	},
	TouchAction = {
		PAN_Y: "pan-y",
		PAN_X: "pan-x",
		NONE: "none"
	},
	fieldOfViewIntrinsics = e => ({
		basis: [numberNode(e[$zoomAdjustedFieldOfView] * Math.PI / 180, "rad")],
		keywords: {
			auto: [null]
		}
	}),
	minFieldOfViewIntrinsics = {
		basis: [degreesToRadians(numberNode(25, "deg"))],
		keywords: {
			auto: [null]
		}
	},
	maxFieldOfViewIntrinsics = e => {
		const t = e[$scene];
		return {
			basis: [degreesToRadians(numberNode(45, "deg"))],
			keywords: {
				auto: [numberNode(t.framedFieldOfView, "deg")]
			}
		}
	},
	cameraOrbitIntrinsics = (() => {
		const e = parseExpressions("0deg 75deg 105%")[0].terms,
			t = normalizeUnit(e[0]),
			i = normalizeUnit(e[1]);
		return e => {
			const n = e[$scene].idealCameraDistance;
			return {
				basis: [t, i, numberNode(n, "m")],
				keywords: {
					auto: [null, null, numberNode(105, "%")]
				}
			}
		}
	})(),
	minCameraOrbitIntrinsics = e => {
		const t = MINIMUM_RADIUS_RATIO * e[$scene].idealCameraDistance;
		return {
			basis: [numberNode(-1 / 0, "rad"), numberNode(Math.PI / 8, "rad"), numberNode(t, "m")],
			keywords: {
				auto: [null, null, null]
			}
		}
	},
	maxCameraOrbitIntrinsics = e => {
		const t = cameraOrbitIntrinsics(e),
			i = new StyleEvaluator([], t).evaluate()[2];
		return {
			basis: [numberNode(1 / 0, "rad"), numberNode(Math.PI - Math.PI / 8, "rad"), numberNode(i, "m")],
			keywords: {
				auto: [null, null, null]
			}
		}
	},
	cameraTargetIntrinsics = e => {
		const t = e[$scene].boundingBox.getCenter(new Vector3);
		return {
			basis: [numberNode(t.x, "m"), numberNode(t.y, "m"), numberNode(t.z, "m")],
			keywords: {
				auto: [null, null, null]
			}
		}
	},
	HALF_PI = Math.PI / 2,
	THIRD_PI = Math.PI / 3,
	QUARTER_PI = HALF_PI / 2,
	TAU = 2 * Math.PI,
	$controls = Symbol("controls"),
	$promptElement = Symbol("promptElement"),
	$promptAnimatedContainer = Symbol("promptAnimatedContainer"),
	$deferInteractionPrompt = Symbol("deferInteractionPrompt"),
	$updateAria = Symbol("updateAria"),
	$updateCameraForRadius = Symbol("updateCameraForRadius"),
	$onBlur = Symbol("onBlur"),
	$onFocus = Symbol("onFocus"),
	$onChange = Symbol("onChange"),
	$onPointerChange = Symbol("onPointerChange"),
	$waitingToPromptUser = Symbol("waitingToPromptUser"),
	$userHasInteracted = Symbol("userHasInteracted"),
	$promptElementVisibleTime = Symbol("promptElementVisibleTime"),
	$lastPromptOffset = Symbol("lastPromptOffset"),
	$focusedTime = Symbol("focusedTime"),
	$zoomAdjustedFieldOfView = Symbol("zoomAdjustedFieldOfView"),
	$lastSpherical = Symbol("lastSpherical"),
	$jumpCamera = Symbol("jumpCamera"),
	$initialized = Symbol("initialized"),
	$maintainThetaPhi = Symbol("maintainThetaPhi"),
	$syncCameraOrbit = Symbol("syncCameraOrbit"),
	$syncFieldOfView = Symbol("syncFieldOfView"),
	$syncCameraTarget = Symbol("syncCameraTarget"),
	$syncMinCameraOrbit = Symbol("syncMinCameraOrbit"),
	$syncMaxCameraOrbit = Symbol("syncMaxCameraOrbit"),
	$syncMinFieldOfView = Symbol("syncMinFieldOfView"),
	$syncMaxFieldOfView = Symbol("syncMaxFieldOfView"),
	ControlsMixin = e => {
		var t, i, n, r, a, s, o, l, c, h, u, d, p, m, A, g, f;
		class v extends e {
			constructor() {
				super(...arguments), this.cameraControls = !1, this.cameraOrbit = "0deg 75deg 105%", this.cameraTarget = "auto auto auto", this.fieldOfView = "auto", this.minCameraOrbit = "auto", this.maxCameraOrbit = "auto", this.minFieldOfView = "auto", this.maxFieldOfView = "auto", this.interactionPromptThreshold = 3e3, this.interactionPromptStyle = InteractionPromptStyle.WIGGLE, this.interactionPrompt = InteractionPromptStrategy.AUTO, this.interactionPolicy = InteractionPolicy.ALWAYS_ALLOW, this.orbitSensitivity = 1, this.touchAction = TouchAction.PAN_Y, this.disableZoom = !1, this.bounds = "legacy", this[t] = this.shadowRoot.querySelector(".interaction-prompt"), this[i] = this.shadowRoot.querySelector(".interaction-prompt > .animated-container"), this[n] = 1 / 0, this[r] = 0, this[a] = 1 / 0, this[s] = !1, this[o] = !1, this[l] = new SmoothControls(this[$scene].camera, this[$userInputElement]), this[c] = 0, this[h] = new Spherical, this[u] = !1, this[d] = !1, this[p] = !1, this[m] = () => {
					const e = this[$userInputElement];
					isFinite(this[$focusedTime]) || (this[$focusedTime] = performance.now());
					const t = this[$ariaLabel];
					e.getAttribute("aria-label") !== t && e.setAttribute("aria-label", t), this.interactionPrompt !== InteractionPromptStrategy.WHEN_FOCUSED || this[$userHasInteracted] || (this[$waitingToPromptUser] = !0)
				}, this[A] = () => {
					this.interactionPrompt === InteractionPromptStrategy.WHEN_FOCUSED && (this[$waitingToPromptUser] = !1, this[$promptElement].classList.remove("visible"), this[$promptElementVisibleTime] = 1 / 0, this[$focusedTime] = 1 / 0)
				}, this[g] = ({
					source: e
				}) => {
					this[$updateAria](), this[$needsRender](), e === ChangeSource.USER_INTERACTION && (this[$userHasInteracted] = !0, this[$deferInteractionPrompt]()), this.dispatchEvent(new CustomEvent("camera-change", {
						detail: {
							source: e
						}
					}))
				}, this[f] = e => {
					"pointer-change-start" === e.type ? this[$container].classList.add("pointer-tumbling") : this[$container].classList.remove("pointer-tumbling")
				}
			}
			getCameraOrbit() {
				const {
					theta: e,
					phi: t,
					radius: i
				} = this[$lastSpherical];
				return {
					theta: e,
					phi: t,
					radius: i,
					toString() {
						return `${this.theta}rad ${this.phi}rad ${this.radius}m`
					}
				}
			}
			getCameraTarget() {
				return toVector3D(this[$renderer].isPresenting ? this[$renderer].arRenderer.target : this[$scene].getTarget())
			}
			getFieldOfView() {
				return this[$controls].getFieldOfView()
			}
			getMinimumFieldOfView() {
				return this[$controls].options.minimumFieldOfView
			}
			getMaximumFieldOfView() {
				return this[$controls].options.maximumFieldOfView
			}
			jumpCameraToGoal() {
				this[$jumpCamera] = !0, this.requestUpdate($jumpCamera, !1)
			}
			resetInteractionPrompt() {
				this[$lastPromptOffset] = 0, this[$promptElementVisibleTime] = 1 / 0, this[$userHasInteracted] = !1, this[$waitingToPromptUser] = this.interactionPrompt === InteractionPromptStrategy.AUTO && this.cameraControls
			}
			connectedCallback() {
				super.connectedCallback(), this[$controls].addEventListener("change", this[$onChange]), this[$controls].addEventListener("pointer-change-start", this[$onPointerChange]), this[$controls].addEventListener("pointer-change-end", this[$onPointerChange])
			}
			disconnectedCallback() {
				super.disconnectedCallback(), this[$controls].removeEventListener("change", this[$onChange]), this[$controls].removeEventListener("pointer-change-start", this[$onPointerChange]), this[$controls].removeEventListener("pointer-change-end", this[$onPointerChange])
			}
			updated(e) {
				super.updated(e);
				const t = this[$controls],
					i = this[$userInputElement];
				if (e.has("cameraControls") && (this.cameraControls ? (t.enableInteraction(), this.interactionPrompt === InteractionPromptStrategy.AUTO && (this[$waitingToPromptUser] = !0), i.addEventListener("focus", this[$onFocus]), i.addEventListener("blur", this[$onBlur])) : (i.removeEventListener("focus", this[$onFocus]), i.removeEventListener("blur", this[$onBlur]), t.disableInteraction(), this[$deferInteractionPrompt]())), e.has("disableZoom") && (t.disableZoom = this.disableZoom), e.has("bounds") && (this[$scene].tightBounds = "tight" === this.bounds), (e.has("interactionPrompt") || e.has("cameraControls") || e.has("src")) && (this.interactionPrompt === InteractionPromptStrategy.AUTO && this.cameraControls && !this[$userHasInteracted] ? this[$waitingToPromptUser] = !0 : this[$deferInteractionPrompt]()), e.has("interactionPromptStyle") && this[$promptElement].classList.toggle("wiggle", this.interactionPromptStyle === InteractionPromptStyle.WIGGLE), e.has("interactionPolicy")) {
					const e = this.interactionPolicy;
					t.applyOptions({
						interactionPolicy: e
					})
				}
				if (e.has("touchAction")) {
					const e = this.touchAction;
					t.applyOptions({
						touchAction: e
					})
				}
				e.has("orbitSensitivity") && (t.sensitivity = this.orbitSensitivity), !0 === this[$jumpCamera] && Promise.resolve().then(() => {
					t.jumpToGoal(), this[$scene].jumpToGoal(), this[$jumpCamera] = !1
				})
			}
			async updateFraming() {
				const e = this[$scene],
					t = e.framedFieldOfView;
				await this.requestUpdate("cameraTarget"), e.updateFraming("tight" === this.bounds ? e.getTarget() : void 0), e.frameModel();
				const i = e.framedFieldOfView,
					n = this[$controls].getFieldOfView() / t;
				this[$zoomAdjustedFieldOfView] = i * n, this[$maintainThetaPhi] = !0, this.requestUpdate("maxFieldOfView"), this.requestUpdate("fieldOfView"), this.requestUpdate("minCameraOrbit"), this.requestUpdate("maxCameraOrbit"), await this.requestUpdate("cameraOrbit")
			} [(t = $promptElement, i = $promptAnimatedContainer, n = $focusedTime, r = $lastPromptOffset, a = $promptElementVisibleTime, s = $userHasInteracted, o = $waitingToPromptUser, l = $controls, c = $zoomAdjustedFieldOfView, h = $lastSpherical, u = $jumpCamera, d = $initialized, p = $maintainThetaPhi, $syncFieldOfView)](e) {
				this[$controls].setFieldOfView(180 * e[0] / Math.PI)
			} [$syncCameraOrbit](e) {
				if (this[$maintainThetaPhi]) {
					const {
						theta: t,
						phi: i
					} = this.getCameraOrbit();
					e[0] = t, e[1] = i, this[$maintainThetaPhi] = !1
				}
				this[$controls].setOrbit(e[0], e[1], e[2])
			} [$syncMinCameraOrbit](e) {
				this[$controls].applyOptions({
					minimumAzimuthalAngle: e[0],
					minimumPolarAngle: e[1],
					minimumRadius: e[2]
				}), this.jumpCameraToGoal()
			} [$syncMaxCameraOrbit](e) {
				this[$controls].applyOptions({
					maximumAzimuthalAngle: e[0],
					maximumPolarAngle: e[1],
					maximumRadius: e[2]
				}), this[$updateCameraForRadius](e[2]), this.jumpCameraToGoal()
			} [$syncMinFieldOfView](e) {
				this[$controls].applyOptions({
					minimumFieldOfView: 180 * e[0] / Math.PI
				}), this.jumpCameraToGoal()
			} [$syncMaxFieldOfView](e) {
				this[$controls].applyOptions({
					maximumFieldOfView: 180 * e[0] / Math.PI
				}), this.jumpCameraToGoal()
			} [$syncCameraTarget](e) {
				const [t, i, n] = e;
				this[$scene].setTarget(t, i, n), this[$renderer].arRenderer.updateTarget()
			} [$tick](e, t) {
				if (super[$tick](e, t), this[$renderer].isPresenting || !this[$hasTransitioned]()) return;
				const i = performance.now();
				if (this[$waitingToPromptUser]) {
					const e = this.interactionPrompt === InteractionPromptStrategy.AUTO ? this[$loadedTime] : this[$focusedTime];
					this.loaded && i > e + this.interactionPromptThreshold && (this[$userInputElement].setAttribute("aria-label", INTERACTION_PROMPT), this[$waitingToPromptUser] = !1, this[$promptElementVisibleTime] = i, this[$promptElement].classList.add("visible"))
				}
				if (isFinite(this[$promptElementVisibleTime]) && this.interactionPromptStyle === InteractionPromptStyle.WIGGLE) {
					const e = this[$scene],
						t = (i - this[$promptElementVisibleTime]) / 5e3 % 1,
						n = wiggle(t),
						r = fade(t);
					if (this[$promptAnimatedContainer].style.opacity = "" + r, n !== this[$lastPromptOffset]) {
						const t = n * e.width * .05,
							i = (n - this[$lastPromptOffset]) * Math.PI / 16;
						this[$promptAnimatedContainer].style.transform = `translateX(${t}px)`, this[$controls].adjustOrbit(i, 0, 0), this[$lastPromptOffset] = n
					}
				}
				this[$controls].update(e, t), this[$scene].updateTarget(t)
			} [$deferInteractionPrompt]() {
				this[$waitingToPromptUser] = !1, this[$promptElement].classList.remove("visible"), this[$promptElementVisibleTime] = 1 / 0
			} [$updateCameraForRadius](e) {
				const {
					idealCameraDistance: t
				} = this[$scene], i = 2 * Math.max(t, e);
				this[$controls].updateNearFar(0, i)
			} [$updateAria]() {
				const {
					theta: e,
					phi: t
				} = this[$lastSpherical], {
					theta: i,
					phi: n
				} = this[$controls].getCameraSpherical(this[$lastSpherical]), r = this.getRootNode();
				if (null != r && r.activeElement === this) {
					const r = (4 + Math.floor((e % TAU + QUARTER_PI) / HALF_PI)) % 4,
						a = (4 + Math.floor((i % TAU + QUARTER_PI) / HALF_PI)) % 4,
						s = Math.floor(t / THIRD_PI),
						o = Math.floor(n / THIRD_PI);
					if (a !== r || o !== s) {
						const e = `View from stage ${POLAR_TRIENT_LABELS[o]}${AZIMUTHAL_QUADRANT_LABELS[a]}`;
						this[$userInputElement].setAttribute("aria-label", e)
					}
				}
			} [$onResize](e) {
				const t = this[$controls],
					i = this[$scene].framedFieldOfView;
				super[$onResize](e);
				const n = this[$scene].framedFieldOfView,
					r = t.getFieldOfView() / i;
				this[$zoomAdjustedFieldOfView] = n * r, t.updateAspect(this[$scene].aspect), this.requestUpdate("maxFieldOfView", this.maxFieldOfView), this.requestUpdate("fieldOfView", this.fieldOfView), this.jumpCameraToGoal()
			} [$onModelLoad]() {
				super[$onModelLoad]();
				const {
					framedFieldOfView: e
				} = this[$scene];
				this[$zoomAdjustedFieldOfView] = e, this[$initialized] ? this[$maintainThetaPhi] = !0 : this[$initialized] = !0, this.requestUpdate("maxFieldOfView", this.maxFieldOfView), this.requestUpdate("fieldOfView", this.fieldOfView), this.requestUpdate("minCameraOrbit", this.minCameraOrbit), this.requestUpdate("maxCameraOrbit", this.maxCameraOrbit), this.requestUpdate("cameraOrbit", this.cameraOrbit), this.requestUpdate("cameraTarget", this.cameraTarget), this.jumpCameraToGoal()
			}
		}
		return m = $onFocus, A = $onBlur, g = $onChange, f = $onPointerChange, __decorate$4([property({
			type: Boolean,
			attribute: "camera-controls"
		})], v.prototype, "cameraControls", void 0), __decorate$4([style({
			intrinsics: cameraOrbitIntrinsics,
			observeEffects: !0,
			updateHandler: $syncCameraOrbit
		}), property({
			type: String,
			attribute: "camera-orbit",
			hasChanged: () => !0
		})], v.prototype, "cameraOrbit", void 0), __decorate$4([style({
			intrinsics: cameraTargetIntrinsics,
			observeEffects: !0,
			updateHandler: $syncCameraTarget
		}), property({
			type: String,
			attribute: "camera-target",
			hasChanged: () => !0
		})], v.prototype, "cameraTarget", void 0), __decorate$4([style({
			intrinsics: fieldOfViewIntrinsics,
			observeEffects: !0,
			updateHandler: $syncFieldOfView
		}), property({
			type: String,
			attribute: "field-of-view",
			hasChanged: () => !0
		})], v.prototype, "fieldOfView", void 0), __decorate$4([style({
			intrinsics: minCameraOrbitIntrinsics,
			updateHandler: $syncMinCameraOrbit
		}), property({
			type: String,
			attribute: "min-camera-orbit",
			hasChanged: () => !0
		})], v.prototype, "minCameraOrbit", void 0), __decorate$4([style({
			intrinsics: maxCameraOrbitIntrinsics,
			updateHandler: $syncMaxCameraOrbit
		}), property({
			type: String,
			attribute: "max-camera-orbit",
			hasChanged: () => !0
		})], v.prototype, "maxCameraOrbit", void 0), __decorate$4([style({
			intrinsics: minFieldOfViewIntrinsics,
			updateHandler: $syncMinFieldOfView
		}), property({
			type: String,
			attribute: "min-field-of-view",
			hasChanged: () => !0
		})], v.prototype, "minFieldOfView", void 0), __decorate$4([style({
			intrinsics: maxFieldOfViewIntrinsics,
			updateHandler: $syncMaxFieldOfView
		}), property({
			type: String,
			attribute: "max-field-of-view",
			hasChanged: () => !0
		})], v.prototype, "maxFieldOfView", void 0), __decorate$4([property({
			type: Number,
			attribute: "interaction-prompt-threshold"
		})], v.prototype, "interactionPromptThreshold", void 0), __decorate$4([property({
			type: String,
			attribute: "interaction-prompt-style"
		})], v.prototype, "interactionPromptStyle", void 0), __decorate$4([property({
			type: String,
			attribute: "interaction-prompt"
		})], v.prototype, "interactionPrompt", void 0), __decorate$4([property({
			type: String,
			attribute: "interaction-policy"
		})], v.prototype, "interactionPolicy", void 0), __decorate$4([property({
			type: Number,
			attribute: "orbit-sensitivity"
		})], v.prototype, "orbitSensitivity", void 0), __decorate$4([property({
			type: String,
			attribute: "touch-action"
		})], v.prototype, "touchAction", void 0), __decorate$4([property({
			type: Boolean,
			attribute: "disable-zoom"
		})], v.prototype, "disableZoom", void 0), __decorate$4([property({
			type: String,
			attribute: "bounds"
		})], v.prototype, "bounds", void 0), v
	};
var __decorate$3 = function(e, t, i, n) {
	for (var r, a = arguments.length, s = a < 3 ? t : null === n ? n = Object.getOwnPropertyDescriptor(t, i) : n, o = e.length - 1; o >= 0; o--)(r = e[o]) && (s = (a < 3 ? r(s) : a > 3 ? r(t, i, s) : r(t, i)) || s);
	return a > 3 && s && Object.defineProperty(t, i, s), s
};
const BASE_OPACITY = .1,
	DEFAULT_SHADOW_INTENSITY = 0,
	DEFAULT_SHADOW_SOFTNESS = 1,
	DEFAULT_EXPOSURE = 1,
	$currentEnvironmentMap = Symbol("currentEnvironmentMap"),
	$applyEnvironmentMap = Symbol("applyEnvironmentMap"),
	$updateEnvironment = Symbol("updateEnvironment"),
	$cancelEnvironmentUpdate = Symbol("cancelEnvironmentUpdate"),
	$onPreload = Symbol("onPreload"),
	EnvironmentMixin = e => {
		var t, i, n;
		class r extends e {
			constructor() {
				super(...arguments), this.environmentImage = null, this.skyboxImage = null, this.shadowIntensity = 0, this.shadowSoftness = 1, this.exposure = 1, this[t] = null, this[i] = null, this[n] = e => {
					e.element === this && this[$updateEnvironment]()
				}
			}
			connectedCallback() {
				super.connectedCallback(), this[$renderer].loader.addEventListener("preload", this[$onPreload])
			}
			disconnectedCallback() {
				super.disconnectedCallback(), this[$renderer].loader.removeEventListener("preload", this[$onPreload])
			}
			updated(e) {
				super.updated(e), e.has("shadowIntensity") && (this[$scene].setShadowIntensity(.1 * this.shadowIntensity), this[$needsRender]()), e.has("shadowSoftness") && (this[$scene].setShadowSoftness(this.shadowSoftness), this[$needsRender]()), e.has("exposure") && (this[$scene].exposure = this.exposure, this[$needsRender]()), (e.has("environmentImage") || e.has("skyboxImage")) && this[$shouldAttemptPreload]() && this[$updateEnvironment]()
			} [(t = $currentEnvironmentMap, i = $cancelEnvironmentUpdate, n = $onPreload, $onModelLoad)]() {
				super[$onModelLoad](), null != this[$currentEnvironmentMap] && this[$applyEnvironmentMap](this[$currentEnvironmentMap])
			}
			async [$updateEnvironment]() {
				const {
					skyboxImage: e,
					environmentImage: t
				} = this;
				null != this[$cancelEnvironmentUpdate] && (this[$cancelEnvironmentUpdate](), this[$cancelEnvironmentUpdate] = null);
				const {
					textureUtils: i
				} = this[$renderer];
				if (null != i) try {
					const {
						environmentMap: n,
						skybox: r
					} = await new Promise(async (n, r) => {
						const a = i.generateEnvironmentMapAndSkybox(deserializeUrl(e), t, {
							progressTracker: this[$progressTracker]
						});
						this[$cancelEnvironmentUpdate] = () => r(a), n(await a)
					}), a = n.texture;
					this[$scene].background = null != r ? r.userData.url === a.userData.url ? a : r : null, this[$applyEnvironmentMap](n.texture), this[$scene].dispatchEvent({
						type: "envmap-update"
					})
				} catch (e) {
					if (e instanceof Error) throw this[$applyEnvironmentMap](null), e
				}
			} [$applyEnvironmentMap](e) {
				this[$currentEnvironmentMap] = e, this[$scene].environment = this[$currentEnvironmentMap], this.dispatchEvent(new CustomEvent("environment-change")), this[$needsRender]()
			}
		}
		return __decorate$3([property({
			type: String,
			attribute: "environment-image"
		})], r.prototype, "environmentImage", void 0), __decorate$3([property({
			type: String,
			attribute: "skybox-image"
		})], r.prototype, "skyboxImage", void 0), __decorate$3([property({
			type: Number,
			attribute: "shadow-intensity"
		})], r.prototype, "shadowIntensity", void 0), __decorate$3([property({
			type: Number,
			attribute: "shadow-softness"
		})], r.prototype, "shadowSoftness", void 0), __decorate$3([property({
			type: Number
		})], r.prototype, "exposure", void 0), r
	};
var _a$4, _b$3;
const INITIAL_STATUS_ANNOUNCEMENT = "This page includes one or more 3D models that are loading",
	FINISHED_LOADING_ANNOUNCEMENT = "All 3D models in the page have loaded",
	UPDATE_STATUS_DEBOUNCE_MS = 100,
	$modelViewerStatusInstance = Symbol("modelViewerStatusInstance"),
	$updateStatus = Symbol("updateStatus");
class LoadingStatusAnnouncer extends EventDispatcher {
	constructor() {
		super(), this[_a$4] = null, this.registeredInstanceStatuses = new Map, this.loadingPromises = [], this.statusElement = document.createElement("p"), this.statusUpdateInProgress = !1, this[_b$3] = debounce(() => this.updateStatus(), 100);
		const {
			statusElement: e
		} = this, {
			style: t
		} = e;
		e.setAttribute("role", "status"), e.classList.add("screen-reader-only"), t.top = t.left = "0", t.pointerEvents = "none"
	}
	registerInstance(e) {
		if (this.registeredInstanceStatuses.has(e)) return;
		let t = () => {};
		const i = !1 === e.loaded && !!e.src,
			n = new Promise(n => {
				if (!i) return void n();
				const r = () => {
					n(), e.removeEventListener("load", r), e.removeEventListener("error", r)
				};
				e.addEventListener("load", r), e.addEventListener("error", r), t = r
			});
		this.registeredInstanceStatuses.set(e, {
			onUnregistered: t
		}), this.loadingPromises.push(n), null == this.modelViewerStatusInstance && (this.modelViewerStatusInstance = e)
	}
	unregisterInstance(e) {
		if (!this.registeredInstanceStatuses.has(e)) return;
		const t = this.registeredInstanceStatuses,
			i = t.get(e);
		t.delete(e), i.onUnregistered(), this.modelViewerStatusInstance === e && (this.modelViewerStatusInstance = t.size > 0 ? getFirstMapKey(t) : null)
	}
	get modelViewerStatusInstance() {
		return this[$modelViewerStatusInstance]
	}
	set modelViewerStatusInstance(e) {
		if (this[$modelViewerStatusInstance] === e) return;
		const {
			statusElement: t
		} = this;
		null != e && null != e.shadowRoot ? e.shadowRoot.appendChild(t) : null != t.parentNode && t.parentNode.removeChild(t), this[$modelViewerStatusInstance] = e, this[$updateStatus]()
	}
	async updateStatus() {
		if (!this.statusUpdateInProgress && 0 !== this.loadingPromises.length) {
			for (this.statusElement.textContent = INITIAL_STATUS_ANNOUNCEMENT, this.statusUpdateInProgress = !0, this.dispatchEvent({
					type: "initial-status-announced"
				}); this.loadingPromises.length;) {
				const {
					loadingPromises: e
				} = this;
				this.loadingPromises = [], await Promise.all(e)
			}
			this.statusElement.textContent = FINISHED_LOADING_ANNOUNCEMENT, this.statusUpdateInProgress = !1, this.dispatchEvent({
				type: "finished-loading-announced"
			})
		}
	}
}
_a$4 = $modelViewerStatusInstance, _b$3 = $updateStatus;
var __decorate$2 = function(e, t, i, n) {
	for (var r, a = arguments.length, s = a < 3 ? t : null === n ? n = Object.getOwnPropertyDescriptor(t, i) : n, o = e.length - 1; o >= 0; o--)(r = e[o]) && (s = (a < 3 ? r(s) : a > 3 ? r(t, i, s) : r(t, i)) || s);
	return a > 3 && s && Object.defineProperty(t, i, s), s
};
const PROGRESS_BAR_UPDATE_THRESHOLD = 100,
	PROGRESS_MASK_BASE_OPACITY = .2,
	DEFAULT_DRACO_DECODER_LOCATION = "https://www.gstatic.com/draco/versioned/decoders/1.3.6/",
	DEFAULT_KTX2_TRANSCODER_LOCATION = "https://www.gstatic.com/basis-universal/versioned/2020-12-21-ef70ddd/",
	SPACE_KEY = 32,
	ENTER_KEY = 13,
	RevealStrategy = {
		AUTO: "auto",
		INTERACTION: "interaction",
		MANUAL: "manual"
	},
	LoadingStrategy = {
		AUTO: "auto",
		LAZY: "lazy",
		EAGER: "eager"
	},
	PosterDismissalSource = {
		INTERACTION: "interaction"
	},
	loadingStatusAnnouncer = new LoadingStatusAnnouncer,
	$defaultProgressBarElement = Symbol("defaultProgressBarElement"),
	$defaultProgressMaskElement = Symbol("defaultProgressMaskElement"),
	$posterContainerElement = Symbol("posterContainerElement"),
	$defaultPosterElement = Symbol("defaultPosterElement"),
	$posterDismissalSource = Symbol("posterDismissalSource"),
	$hidePoster = Symbol("hidePoster"),
	$modelIsRevealed = Symbol("modelIsRevealed"),
	$updateProgressBar = Symbol("updateProgressBar"),
	$lastReportedProgress = Symbol("lastReportedProgress"),
	$transitioned = Symbol("transitioned"),
	$ariaLabelCallToAction = Symbol("ariaLabelCallToAction"),
	$onClick = Symbol("onClick"),
	$onKeydown = Symbol("onKeydown"),
	$onProgress = Symbol("onProgress"),
	LoadingMixin = e => {
		var t, i, n, r, a, s, o, l, c, h, u, d, p;
		class m extends e {
			constructor(...e) {
				super(...e), this.poster = null, this.reveal = RevealStrategy.AUTO, this.loading = LoadingStrategy.AUTO, this[t] = !1, this[i] = !1, this[n] = 0, this[r] = null, this[a] = this.shadowRoot.querySelector(".slot.poster"), this[s] = this.shadowRoot.querySelector("#default-poster"), this[o] = this.shadowRoot.querySelector("#default-progress-bar > .bar"), this[l] = this.shadowRoot.querySelector("#default-progress-bar > .mask"), this[c] = this[$defaultPosterElement].getAttribute("aria-label"), this[h] = throttle(e => {
					const t = this[$defaultProgressBarElement].parentNode;
					requestAnimationFrame(() => {
						this[$defaultProgressMaskElement].style.opacity = "" + .2 * (1 - e), this[$defaultProgressBarElement].style.transform = `scaleX(${e})`, 0 === e && (t.removeChild(this[$defaultProgressBarElement]), t.appendChild(this[$defaultProgressBarElement])), 1 === e ? this[$defaultProgressBarElement].classList.add("hide") : this[$defaultProgressBarElement].classList.remove("hide")
					})
				}, 100), this[u] = () => {
					this.reveal !== RevealStrategy.MANUAL && this.dismissPoster()
				}, this[d] = e => {
					if (this.reveal !== RevealStrategy.MANUAL) switch (e.keyCode) {
						case 32:
						case 13:
							this.dismissPoster()
					}
				}, this[p] = e => {
					const t = e.detail.totalProgress;
					this[$lastReportedProgress] = Math.max(t, this[$lastReportedProgress]), 1 === t && (this[$updateProgressBar].flush(), !this[$sceneIsReady]() || null == this[$posterDismissalSource] && this.reveal !== RevealStrategy.AUTO || this[$hidePoster]()), this[$updateProgressBar](t), this.dispatchEvent(new CustomEvent("progress", {
						detail: {
							totalProgress: t
						}
					}))
				};
				const m = self.ModelViewerElement || {},
					A = m.dracoDecoderLocation || DEFAULT_DRACO_DECODER_LOCATION;
				CachingGLTFLoader.setDRACODecoderLocation(A);
				const g = m.ktx2TranscoderLocation || DEFAULT_KTX2_TRANSCODER_LOCATION;
				CachingGLTFLoader.setKTX2TranscoderLocation(g)
			}
			static set dracoDecoderLocation(e) {
				CachingGLTFLoader.setDRACODecoderLocation(e)
			}
			static get dracoDecoderLocation() {
				return CachingGLTFLoader.getDRACODecoderLocation()
			}
			static set ktx2TranscoderLocation(e) {
				CachingGLTFLoader.setKTX2TranscoderLocation(e)
			}
			static get ktx2TranscoderLocation() {
				return CachingGLTFLoader.getKTX2TranscoderLocation()
			}
			static mapURLs(e) {
				Renderer.singleton.loader[$loader].manager.setURLModifier(e)
			}
			dismissPoster() {
				this[$sceneIsReady]() ? this[$hidePoster]() : (this[$posterDismissalSource] = PosterDismissalSource.INTERACTION, this[$updateSource]())
			}
			showPoster() {
				const e = this[$posterContainerElement],
					t = this[$defaultPosterElement];
				t.removeAttribute("tabindex"), t.removeAttribute("aria-hidden"), e.classList.add("show");
				const i = this.modelIsVisible;
				this[$modelIsRevealed] = !1, this[$announceModelVisibility](i), this[$transitioned] = !1
			}
			getDimensions() {
				return toVector3D(this[$scene].size)
			}
			connectedCallback() {
				super.connectedCallback(), this[$posterContainerElement].addEventListener("click", this[$onClick]), this[$posterContainerElement].addEventListener("keydown", this[$onKeydown]), this[$progressTracker].addEventListener("progress", this[$onProgress]), loadingStatusAnnouncer.registerInstance(this)
			}
			disconnectedCallback() {
				super.disconnectedCallback(), this[$posterContainerElement].removeEventListener("click", this[$onClick]), this[$posterContainerElement].removeEventListener("keydown", this[$onKeydown]), this[$progressTracker].removeEventListener("progress", this[$onProgress]), loadingStatusAnnouncer.unregisterInstance(this)
			}
			async updated(e) {
				super.updated(e), e.has("poster") && null != this.poster && (this[$defaultPosterElement].style.backgroundImage = `url(${this.poster})`), e.has("alt") && this[$defaultPosterElement].setAttribute("aria-label", `${this[$ariaLabel]}. ${this[$ariaLabelCallToAction]}`), (e.has("reveal") || e.has("loaded")) && (this[$sceneIsReady]() || this[$updateSource]())
			} [(t = $modelIsRevealed, i = $transitioned, n = $lastReportedProgress, r = $posterDismissalSource, a = $posterContainerElement, s = $defaultPosterElement, o = $defaultProgressBarElement, l = $defaultProgressMaskElement, c = $ariaLabelCallToAction, h = $updateProgressBar, u = $onClick, d = $onKeydown, p = $onProgress, $shouldAttemptPreload)]() {
				return !!this.src && (null != this[$posterDismissalSource] || this.loading === LoadingStrategy.EAGER || this.reveal === RevealStrategy.AUTO && this[$isElementInViewport])
			} [$sceneIsReady]() {
				const {
					src: e
				} = this;
				return !!e && super[$sceneIsReady]() && 1 === this[$lastReportedProgress]
			} [$hidePoster]() {
				this[$posterDismissalSource] = null;
				const e = this[$posterContainerElement],
					t = this[$defaultPosterElement];
				if (e.classList.contains("show")) {
					e.classList.remove("show");
					const i = this.modelIsVisible;
					this[$modelIsRevealed] = !0, this[$announceModelVisibility](i), e.addEventListener("transitionend", () => {
						requestAnimationFrame(() => {
							this[$transitioned] = !0;
							const e = this.getRootNode();
							e && e.activeElement === this && this[$userInputElement].focus(), t.setAttribute("aria-hidden", "true"), t.tabIndex = -1, this.dispatchEvent(new CustomEvent("poster-dismissed"))
						})
					}, {
						once: !0
					})
				}
			} [$getModelIsVisible]() {
				return super[$getModelIsVisible]() && this[$modelIsRevealed]
			} [$hasTransitioned]() {
				return super[$hasTransitioned]() && this[$transitioned]
			}
			async [$updateSource]() {
				this[$lastReportedProgress] = 0, null != this[$scene].currentGLTF && null != this.src && this[$shouldAttemptPreload]() || this.showPoster(), await super[$updateSource]()
			}
		}
		return __decorate$2([property({
			type: String
		})], m.prototype, "poster", void 0), __decorate$2([property({
			type: String
		})], m.prototype, "reveal", void 0), __decorate$2([property({
			type: String
		})], m.prototype, "loading", void 0), m
	};
var GLTFExporter = function() {
	function e() {
		this.pluginCallbacks = [], this.register((function(e) {
			return new M(e)
		})), this.register((function(e) {
			return new T(e)
		})), this.register((function(e) {
			return new B(e)
		}))
	}
	e.prototype = {
		constructor: e,
		register: function(e) {
			return -1 === this.pluginCallbacks.indexOf(e) && this.pluginCallbacks.push(e), this
		},
		unregister: function(e) {
			return -1 !== this.pluginCallbacks.indexOf(e) && this.pluginCallbacks.splice(this.pluginCallbacks.indexOf(e), 1), this
		},
		parse: function(e, t, i) {
			for (var n = new I, r = [], a = 0, s = this.pluginCallbacks.length; a < s; a++) r.push(this.pluginCallbacks[a](n));
			n.setPlugins(r), n.write(e, t, i)
		}
	};
	var t = 0,
		i = 1,
		n = 2,
		r = 3,
		a = 4,
		s = 5121,
		o = 5123,
		l = 5126,
		c = 5125,
		h = 34962,
		u = 34963,
		d = 9728,
		p = 9729,
		m = 9984,
		A = 9985,
		g = 9986,
		f = 9987,
		v = 33071,
		y = 33648,
		E = 10497,
		_ = {};
	_[1003] = d, _[1004] = m, _[1005] = g, _[1006] = p, _[1007] = A, _[1008] = f, _[1001] = v, _[1e3] = E, _[1002] = y;
	var b = {
		scale: "scale",
		position: "translation",
		quaternion: "rotation",
		morphTargetInfluences: "weights"
	};

	function x(e, t) {
		return e.length === t.length && e.every((function(e, i) {
			return e === t[i]
		}))
	}

	function w(e) {
		return 4 * Math.ceil(e / 4)
	}

	function C(e, t) {
		t = t || 0;
		var i = w(e.byteLength);
		if (i !== e.byteLength) {
			var n = new Uint8Array(i);
			if (n.set(new Uint8Array(e)), 0 !== t)
				for (var r = e.byteLength; r < i; r++) n[r] = t;
			return n.buffer
		}
		return e
	}
	var S = null;

	function I() {
		this.plugins = [], this.options = {}, this.pending = [], this.buffers = [], this.byteOffset = 0, this.buffers = [], this.nodeMap = new Map, this.skins = [], this.extensionsUsed = {}, this.uids = new Map, this.uid = 0, this.json = {
			asset: {
				version: "2.0",
				generator: "THREE.GLTFExporter"
			}
		}, this.cache = {
			meshes: new Map,
			attributes: new Map,
			attributesNormalized: new Map,
			materials: new Map,
			textures: new Map,
			images: new Map
		}
	}

	function M(e) {
		this.writer = e, this.name = "KHR_lights_punctual"
	}

	function T(e) {
		this.writer = e, this.name = "KHR_materials_unlit"
	}

	function B(e) {
		this.writer = e, this.name = "KHR_materials_pbrSpecularGlossiness"
	}
	return I.prototype = {
		constructor: I,
		setPlugins: function(e) {
			this.plugins = e
		},
		write: function(e, t, i) {
			this.options = Object.assign({}, {
				binary: !1,
				trs: !1,
				onlyVisible: !0,
				truncateDrawRange: !0,
				embedImages: !0,
				maxTextureSize: 1 / 0,
				animations: [],
				includeCustomExtensions: !1
			}, i), this.options.animations.length > 0 && (this.options.trs = !0), this.processInput(e);
			var n = this;
			Promise.all(this.pending).then((function() {
				var e, i = n.buffers,
					r = n.json,
					a = n.options,
					s = n.extensionsUsed,
					o = new Blob(i, {
						type: "application/octet-stream"
					}),
					l = Object.keys(s);
				(l.length > 0 && (r.extensionsUsed = l), r.buffers && r.buffers.length > 0 && (r.buffers[0].byteLength = o.size), !0 === a.binary) ? ((e = new window.FileReader).readAsArrayBuffer(o), e.onloadend = function() {
					var i = C(e.result),
						n = new DataView(new ArrayBuffer(8));
					n.setUint32(0, i.byteLength, !0), n.setUint32(4, 5130562, !0);
					var a = C(function(e) {
							if (void 0 !== window.TextEncoder) return (new TextEncoder).encode(e).buffer;
							for (var t = new Uint8Array(new ArrayBuffer(e.length)), i = 0, n = e.length; i < n; i++) {
								var r = e.charCodeAt(i);
								t[i] = r > 255 ? 32 : r
							}
							return t.buffer
						}(JSON.stringify(r)), 32),
						s = new DataView(new ArrayBuffer(8));
					s.setUint32(0, a.byteLength, !0), s.setUint32(4, 1313821514, !0);
					var o = new ArrayBuffer(12),
						l = new DataView(o);
					l.setUint32(0, 1179937895, !0), l.setUint32(4, 2, !0);
					var c = 12 + s.byteLength + a.byteLength + n.byteLength + i.byteLength;
					l.setUint32(8, c, !0);
					var h = new Blob([o, s, a, n, i], {
							type: "application/octet-stream"
						}),
						u = new window.FileReader;
					u.readAsArrayBuffer(h), u.onloadend = function() {
						t(u.result)
					}
				}) : r.buffers && r.buffers.length > 0 ? ((e = new window.FileReader).readAsDataURL(o), e.onloadend = function() {
					var i = e.result;
					r.buffers[0].uri = i, t(r)
				}) : t(r)
			}))
		},
		serializeUserData: function(e, t) {
			if (0 !== Object.keys(e.userData).length) {
				var i = this.options,
					n = this.extensionsUsed;
				try {
					var r = JSON.parse(JSON.stringify(e.userData));
					if (i.includeCustomExtensions && r.gltfExtensions) {
						for (var a in void 0 === t.extensions && (t.extensions = {}), r.gltfExtensions) t.extensions[a] = r.gltfExtensions[a], n[a] = !0;
						delete r.gltfExtensions
					}
					Object.keys(r).length > 0 && (t.extras = r)
				} catch (t) {
					console.warn("THREE.GLTFExporter: userData of '" + e.name + "' won't be serialized because of JSON.stringify error - " + t.message)
				}
			}
		},
		getUID: function(e) {
			return this.uids.has(e) || this.uids.set(e, this.uid++), this.uids.get(e)
		},
		isNormalizedNormalAttribute: function(e) {
			if (this.cache.attributesNormalized.has(e)) return !1;
			for (var t = new Vector3, i = 0, n = e.count; i < n; i++)
				if (Math.abs(t.fromBufferAttribute(e, i).length() - 1) > 5e-4) return !1;
			return !0
		},
		createNormalizedNormalAttribute: function(e) {
			var t = this.cache;
			if (t.attributesNormalized.has(e)) return t.attributesNormalized.get(e);
			for (var i = e.clone(), n = new Vector3, r = 0, a = i.count; r < a; r++) n.fromBufferAttribute(i, r), 0 === n.x && 0 === n.y && 0 === n.z ? n.setX(1) : n.normalize(), i.setXYZ(r, n.x, n.y, n.z);
			return t.attributesNormalized.set(e, i), i
		},
		applyTextureTransform: function(e, t) {
			var i = !1,
				n = {};
			0 === t.offset.x && 0 === t.offset.y || (n.offset = t.offset.toArray(), i = !0), 0 !== t.rotation && (n.rotation = t.rotation, i = !0), 1 === t.repeat.x && 1 === t.repeat.y || (n.scale = t.repeat.toArray(), i = !0), i && (e.extensions = e.extensions || {}, e.extensions.KHR_texture_transform = n, this.extensionsUsed.KHR_texture_transform = !0)
		},
		processBuffer: function(e) {
			var t = this.json,
				i = this.buffers;
			return t.buffers || (t.buffers = [{
				byteLength: 0
			}]), i.push(e), 0
		},
		processBufferView: function(e, t, i, n, r) {
			var a, u = this.json;
			u.bufferViews || (u.bufferViews = []), a = t === s ? 1 : t === o ? 2 : 4;
			for (var d = w(n * e.itemSize * a), p = new DataView(new ArrayBuffer(d)), m = 0, A = i; A < i + n; A++)
				for (var g = 0; g < e.itemSize; g++) {
					var f;
					e.itemSize > 4 ? f = e.array[A * e.itemSize + g] : 0 === g ? f = e.getX(A) : 1 === g ? f = e.getY(A) : 2 === g ? f = e.getZ(A) : 3 === g && (f = e.getW(A)), t === l ? p.setFloat32(m, f, !0) : t === c ? p.setUint32(m, f, !0) : t === o ? p.setUint16(m, f, !0) : t === s && p.setUint8(m, f), m += a
				}
			var v = {
				buffer: this.processBuffer(p.buffer),
				byteOffset: this.byteOffset,
				byteLength: d
			};
			return void 0 !== r && (v.target = r), r === h && (v.byteStride = e.itemSize * a), this.byteOffset += d, u.bufferViews.push(v), {
				id: u.bufferViews.length - 1,
				byteLength: 0
			}
		},
		processBufferViewImage: function(e) {
			var t = this,
				i = t.json;
			return i.bufferViews || (i.bufferViews = []), new Promise((function(n) {
				var r = new window.FileReader;
				r.readAsArrayBuffer(e), r.onloadend = function() {
					var e = C(r.result),
						a = {
							buffer: t.processBuffer(e),
							byteOffset: t.byteOffset,
							byteLength: e.byteLength
						};
					t.byteOffset += e.byteLength, n(i.bufferViews.push(a) - 1)
				}
			}))
		},
		processAccessor: function(e, t, i, n) {
			var r, a = this.options,
				d = this.json;
			if (e.array.constructor === Float32Array) r = l;
			else if (e.array.constructor === Uint32Array) r = c;
			else if (e.array.constructor === Uint16Array) r = o;
			else {
				if (e.array.constructor !== Uint8Array) throw new Error("THREE.GLTFExporter: Unsupported bufferAttribute component type.");
				r = s
			}
			if (void 0 === i && (i = 0), void 0 === n && (n = e.count), a.truncateDrawRange && void 0 !== t && null === t.index) {
				var p = i + n,
					m = t.drawRange.count === 1 / 0 ? e.count : t.drawRange.start + t.drawRange.count;
				i = Math.max(i, t.drawRange.start), (n = Math.min(p, m) - i) < 0 && (n = 0)
			}
			if (0 === n) return null;
			var A, g = function(e, t, i) {
				for (var n = {
						min: new Array(e.itemSize).fill(Number.POSITIVE_INFINITY),
						max: new Array(e.itemSize).fill(Number.NEGATIVE_INFINITY)
					}, r = t; r < t + i; r++)
					for (var a = 0; a < e.itemSize; a++) {
						var s;
						e.itemSize > 4 ? s = e.array[r * e.itemSize + a] : 0 === a ? s = e.getX(r) : 1 === a ? s = e.getY(r) : 2 === a ? s = e.getZ(r) : 3 === a && (s = e.getW(r)), n.min[a] = Math.min(n.min[a], s), n.max[a] = Math.max(n.max[a], s)
					}
				return n
			}(e, i, n);
			void 0 !== t && (A = e === t.index ? u : h);
			var f = this.processBufferView(e, r, i, n, A),
				v = {
					bufferView: f.id,
					byteOffset: f.byteOffset,
					componentType: r,
					count: n,
					max: g.max,
					min: g.min,
					type: {
						1: "SCALAR",
						2: "VEC2",
						3: "VEC3",
						4: "VEC4",
						16: "MAT4"
					} [e.itemSize]
				};
			return !0 === e.normalized && (v.normalized = !0), d.accessors || (d.accessors = []), d.accessors.push(v) - 1
		},
		processImage: function(e, t, i) {
			var n = this,
				r = n.cache,
				a = n.json,
				s = n.options,
				o = n.pending;
			r.images.has(e) || r.images.set(e, {});
			var l = r.images.get(e),
				c = 1023 === t ? "image/png" : "image/jpeg",
				h = c + ":flipY/" + i.toString();
			if (void 0 !== l[h]) return l[h];
			a.images || (a.images = []);
			var u = {
				mimeType: c
			};
			if (s.embedImages) {
				var d = S = S || document.createElement("canvas");
				d.width = Math.min(e.width, s.maxTextureSize), d.height = Math.min(e.height, s.maxTextureSize);
				var p = d.getContext("2d");
				if (!0 === i && (p.translate(0, d.height), p.scale(1, -1)), "undefined" != typeof HTMLImageElement && e instanceof HTMLImageElement || "undefined" != typeof HTMLCanvasElement && e instanceof HTMLCanvasElement || "undefined" != typeof OffscreenCanvas && e instanceof OffscreenCanvas || "undefined" != typeof ImageBitmap && e instanceof ImageBitmap) p.drawImage(e, 0, 0, d.width, d.height);
				else {
					1023 !== t && 1022 !== t && console.error("GLTFExporter: Only RGB and RGBA formats are supported."), (e.width > s.maxTextureSize || e.height > s.maxTextureSize) && console.warn("GLTFExporter: Image size is bigger than maxTextureSize", e);
					var m = e.data;
					if (1022 === t) {
						m = new Uint8ClampedArray(e.height * e.width * 4);
						for (var A = 0, g = 0; A < m.length; A += 4, g += 3) m[A + 0] = e.data[g + 0], m[A + 1] = e.data[g + 1], m[A + 2] = e.data[g + 2], m[A + 3] = 255
					}
					p.putImageData(new ImageData(m, e.width, e.height), 0, 0)
				}!0 === s.binary ? o.push(new Promise((function(e) {
					d.toBlob((function(t) {
						n.processBufferViewImage(t).then((function(t) {
							u.bufferView = t, e()
						}))
					}), c)
				}))) : u.uri = d.toDataURL(c)
			} else u.uri = e.src;
			var f = a.images.push(u) - 1;
			return l[h] = f, f
		},
		processSampler: function(e) {
			var t = this.json;
			t.samplers || (t.samplers = []);
			var i = {
				magFilter: _[e.magFilter],
				minFilter: _[e.minFilter],
				wrapS: _[e.wrapS],
				wrapT: _[e.wrapT]
			};
			return t.samplers.push(i) - 1
		},
		processTexture: function(e) {
			var t = this.cache,
				i = this.json;
			if (t.textures.has(e)) return t.textures.get(e);
			i.textures || (i.textures = []);
			var n = {
				sampler: this.processSampler(e),
				source: this.processImage(e.image, e.format, e.flipY)
			};
			e.name && (n.name = e.name), this._invokeAll((function(t) {
				t.writeTexture && t.writeTexture(e, n)
			}));
			var r = i.textures.push(n) - 1;
			return t.textures.set(e, r), r
		},
		processMaterial: function(e) {
			var t = this.cache,
				i = this.json;
			if (t.materials.has(e)) return t.materials.get(e);
			if (e.isShaderMaterial) return console.warn("GLTFExporter: THREE.ShaderMaterial not supported."), null;
			i.materials || (i.materials = []);
			var n = {
				pbrMetallicRoughness: {}
			};
			!0 !== e.isMeshStandardMaterial && !0 !== e.isMeshBasicMaterial && console.warn("GLTFExporter: Use MeshStandardMaterial or MeshBasicMaterial for best results.");
			var r = e.color.toArray().concat([e.opacity]);
			if (x(r, [1, 1, 1, 1]) || (n.pbrMetallicRoughness.baseColorFactor = r), e.isMeshStandardMaterial ? (n.pbrMetallicRoughness.metallicFactor = e.metalness, n.pbrMetallicRoughness.roughnessFactor = e.roughness) : (n.pbrMetallicRoughness.metallicFactor = .5, n.pbrMetallicRoughness.roughnessFactor = .5), e.metalnessMap || e.roughnessMap)
				if (e.metalnessMap === e.roughnessMap) {
					var a = {
						index: this.processTexture(e.metalnessMap)
					};
					this.applyTextureTransform(a, e.metalnessMap), n.pbrMetallicRoughness.metallicRoughnessTexture = a
				} else console.warn("THREE.GLTFExporter: Ignoring metalnessMap and roughnessMap because they are not the same Texture.");
			if (e.map) {
				var s = {
					index: this.processTexture(e.map)
				};
				this.applyTextureTransform(s, e.map), n.pbrMetallicRoughness.baseColorTexture = s
			}
			if (e.emissive) {
				var o = e.emissive.clone().multiplyScalar(e.emissiveIntensity).toArray();
				if (x(o, [0, 0, 0]) || (n.emissiveFactor = o), e.emissiveMap) {
					var l = {
						index: this.processTexture(e.emissiveMap)
					};
					this.applyTextureTransform(l, e.emissiveMap), n.emissiveTexture = l
				}
			}
			if (e.normalMap) {
				var c = {
					index: this.processTexture(e.normalMap)
				};
				e.normalScale && -1 !== e.normalScale.x && (e.normalScale.x !== e.normalScale.y && console.warn("THREE.GLTFExporter: Normal scale components are different, ignoring Y and exporting X."), c.scale = e.normalScale.x), this.applyTextureTransform(c, e.normalMap), n.normalTexture = c
			}
			if (e.aoMap) {
				var h = {
					index: this.processTexture(e.aoMap),
					texCoord: 1
				};
				1 !== e.aoMapIntensity && (h.strength = e.aoMapIntensity), this.applyTextureTransform(h, e.aoMap), n.occlusionTexture = h
			}
			e.transparent ? n.alphaMode = "BLEND" : e.alphaTest > 0 && (n.alphaMode = "MASK", n.alphaCutoff = e.alphaTest), 2 === e.side && (n.doubleSided = !0), "" !== e.name && (n.name = e.name), this.serializeUserData(e, n), this._invokeAll((function(t) {
				t.writeMaterial && t.writeMaterial(e, n)
			}));
			var u = i.materials.push(n) - 1;
			return t.materials.set(e, u), u
		},
		processMesh: function(e) {
			var s = this.cache,
				o = this.json,
				l = [e.geometry.uuid];
			if (Array.isArray(e.material))
				for (var c = 0, h = e.material.length; c < h; c++) l.push(e.material[c].uuid);
			else l.push(e.material.uuid);
			var u = l.join(":");
			if (s.meshes.has(u)) return s.meshes.get(u);
			var d, p = e.geometry;
			if (d = e.isLineSegments ? i : e.isLineLoop ? n : e.isLine ? r : e.isPoints ? t : e.material.wireframe ? i : a, !0 !== p.isBufferGeometry) throw new Error("THREE.GLTFExporter: Geometry is not of type THREE.BufferGeometry.");
			var m = {},
				A = {},
				g = [],
				f = [],
				v = {
					uv: "TEXCOORD_0",
					uv2: "TEXCOORD_1",
					color: "COLOR_0",
					skinWeight: "WEIGHTS_0",
					skinIndex: "JOINTS_0"
				},
				y = p.getAttribute("normal");
			void 0 === y || this.isNormalizedNormalAttribute(y) || (console.warn("THREE.GLTFExporter: Creating normalized normal attribute from the non-normalized one."), p.setAttribute("normal", this.createNormalizedNormalAttribute(y)));
			var E = null;
			for (var _ in p.attributes)
				if ("morph" !== _.substr(0, 5)) {
					var b = p.attributes[_];
					_ = v[_] || _.toUpperCase();
					if (/^(POSITION|NORMAL|TANGENT|TEXCOORD_\d+|COLOR_\d+|JOINTS_\d+|WEIGHTS_\d+)$/.test(_) || (_ = "_" + _), s.attributes.has(this.getUID(b))) A[_] = s.attributes.get(this.getUID(b));
					else {
						E = null;
						var x = b.array;
						"JOINTS_0" !== _ || x instanceof Uint16Array || x instanceof Uint8Array || (console.warn('GLTFExporter: Attribute "skinIndex" converted to type UNSIGNED_SHORT.'), E = new BufferAttribute(new Uint16Array(x), b.itemSize, b.normalized));
						var w = this.processAccessor(E || b, p);
						null !== w && (A[_] = w, s.attributes.set(this.getUID(b), w))
					}
				} if (void 0 !== y && p.setAttribute("normal", y), 0 === Object.keys(A).length) return null;
			if (void 0 !== e.morphTargetInfluences && e.morphTargetInfluences.length > 0) {
				var C = [],
					S = [],
					I = {};
				if (void 0 !== e.morphTargetDictionary)
					for (var M in e.morphTargetDictionary) I[e.morphTargetDictionary[M]] = M;
				for (c = 0; c < e.morphTargetInfluences.length; ++c) {
					var T = {},
						B = !1;
					for (var _ in p.morphAttributes)
						if ("position" === _ || "normal" === _) {
							b = p.morphAttributes[_][c];
							var L = _.toUpperCase(),
								R = p.attributes[_];
							if (s.attributes.has(this.getUID(b))) T[L] = s.attributes.get(this.getUID(b));
							else {
								var D = b.clone();
								if (!p.morphTargetsRelative)
									for (var P = 0, Q = b.count; P < Q; P++) D.setXYZ(P, b.getX(P) - R.getX(P), b.getY(P) - R.getY(P), b.getZ(P) - R.getZ(P));
								T[L] = this.processAccessor(D, p), s.attributes.set(this.getUID(R), T[L])
							}
						} else B || (console.warn("GLTFExporter: Only POSITION and NORMAL morph are supported."), B = !0);
					f.push(T), C.push(e.morphTargetInfluences[c]), void 0 !== e.morphTargetDictionary && S.push(I[c])
				}
				m.weights = C, S.length > 0 && (m.extras = {}, m.extras.targetNames = S)
			}
			var F = Array.isArray(e.material);
			if (F && 0 === p.groups.length) return null;
			for (var O = F ? e.material : [e.material], N = F ? p.groups : [{
					materialIndex: 0,
					start: void 0,
					count: void 0
				}], k = (c = 0, N.length); c < k; c++) {
				var U = {
					mode: d,
					attributes: A
				};
				if (this.serializeUserData(p, U), f.length > 0 && (U.targets = f), null !== p.index) {
					var G = this.getUID(p.index);
					void 0 === N[c].start && void 0 === N[c].count || (G += ":" + N[c].start + ":" + N[c].count), s.attributes.has(G) ? U.indices = s.attributes.get(G) : (U.indices = this.processAccessor(p.index, p, N[c].start, N[c].count), s.attributes.set(G, U.indices)), null === U.indices && delete U.indices
				}
				var H = this.processMaterial(O[N[c].materialIndex]);
				null !== H && (U.material = H), g.push(U)
			}
			m.primitives = g, o.meshes || (o.meshes = []), this._invokeAll((function(t) {
				t.writeMesh && t.writeMesh(e, m)
			}));
			var V = o.meshes.push(m) - 1;
			return s.meshes.set(u, V), V
		},
		processCamera: function(e) {
			var t = this.json;
			t.cameras || (t.cameras = []);
			var i = e.isOrthographicCamera,
				n = {
					type: i ? "orthographic" : "perspective"
				};
			return i ? n.orthographic = {
				xmag: 2 * e.right,
				ymag: 2 * e.top,
				zfar: e.far <= 0 ? .001 : e.far,
				znear: e.near < 0 ? 0 : e.near
			} : n.perspective = {
				aspectRatio: e.aspect,
				yfov: MathUtils.degToRad(e.fov),
				zfar: e.far <= 0 ? .001 : e.far,
				znear: e.near < 0 ? 0 : e.near
			}, "" !== e.name && (n.name = e.type), t.cameras.push(n) - 1
		},
		processAnimation: function(t, i) {
			var n = this.json,
				r = this.nodeMap;
			n.animations || (n.animations = []);
			for (var a = (t = e.Utils.mergeMorphTargetTracks(t.clone(), i)).tracks, s = [], o = [], l = 0; l < a.length; ++l) {
				var c = a[l],
					h = PropertyBinding.parseTrackName(c.name),
					u = PropertyBinding.findNode(i, h.nodeName),
					d = b[h.propertyName];
				if ("bones" === h.objectName && (u = !0 === u.isSkinnedMesh ? u.skeleton.getBoneByName(h.objectIndex) : void 0), !u || !d) return console.warn('THREE.GLTFExporter: Could not export animation track "%s".', c.name), null;
				var p, m = c.values.length / c.times.length;
				d === b.morphTargetInfluences && (m /= u.morphTargetInfluences.length), !0 === c.createInterpolant.isInterpolantFactoryMethodGLTFCubicSpline ? (p = "CUBICSPLINE", m /= 3) : p = 2300 === c.getInterpolation() ? "STEP" : "LINEAR", o.push({
					input: this.processAccessor(new BufferAttribute(c.times, 1)),
					output: this.processAccessor(new BufferAttribute(c.values, m)),
					interpolation: p
				}), s.push({
					sampler: o.length - 1,
					target: {
						node: r.get(u),
						path: d
					}
				})
			}
			return n.animations.push({
				name: t.name || "clip_" + n.animations.length,
				samplers: o,
				channels: s
			}), n.animations.length - 1
		},
		processSkin: function(e) {
			var t = this.json,
				i = this.nodeMap,
				n = t.nodes[i.get(e)],
				r = e.skeleton;
			if (void 0 === r) return null;
			var a = e.skeleton.bones[0];
			if (void 0 === a) return null;
			for (var s = [], o = new Float32Array(16 * r.bones.length), l = new Matrix4, c = 0; c < r.bones.length; ++c) s.push(i.get(r.bones[c])), l.copy(r.boneInverses[c]), l.multiply(e.bindMatrix).toArray(o, 16 * c);
			return void 0 === t.skins && (t.skins = []), t.skins.push({
				inverseBindMatrices: this.processAccessor(new BufferAttribute(o, 16)),
				joints: s,
				skeleton: i.get(a)
			}), n.skin = t.skins.length - 1
		},
		processNode: function(e) {
			var t = this.json,
				i = this.options,
				n = this.nodeMap;
			t.nodes || (t.nodes = []);
			var r = {};
			if (i.trs) {
				var a = e.quaternion.toArray(),
					s = e.position.toArray(),
					o = e.scale.toArray();
				x(a, [0, 0, 0, 1]) || (r.rotation = a), x(s, [0, 0, 0]) || (r.translation = s), x(o, [1, 1, 1]) || (r.scale = o)
			} else e.matrixAutoUpdate && e.updateMatrix(), !1 === x(e.matrix.elements, [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]) && (r.matrix = e.matrix.elements);
			if ("" !== e.name && (r.name = String(e.name)), this.serializeUserData(e, r), e.isMesh || e.isLine || e.isPoints) {
				var l = this.processMesh(e);
				null !== l && (r.mesh = l)
			} else e.isCamera && (r.camera = this.processCamera(e));
			if (e.isSkinnedMesh && this.skins.push(e), e.children.length > 0) {
				for (var c = [], h = 0, u = e.children.length; h < u; h++) {
					var d = e.children[h];
					if (d.visible || !1 === i.onlyVisible) null !== (p = this.processNode(d)) && c.push(p)
				}
				c.length > 0 && (r.children = c)
			}
			this._invokeAll((function(t) {
				t.writeNode && t.writeNode(e, r)
			}));
			var p = t.nodes.push(r) - 1;
			return n.set(e, p), p
		},
		processScene: function(e) {
			var t = this.json,
				i = this.options;
			t.scenes || (t.scenes = [], t.scene = 0);
			var n = {};
			"" !== e.name && (n.name = e.name), t.scenes.push(n);
			for (var r = [], a = 0, s = e.children.length; a < s; a++) {
				var o = e.children[a];
				if (o.visible || !1 === i.onlyVisible) {
					var l = this.processNode(o);
					null !== l && r.push(l)
				}
			}
			r.length > 0 && (n.nodes = r), this.serializeUserData(e, n)
		},
		processObjects: function(e) {
			var t = new Scene;
			t.name = "AuxScene";
			for (var i = 0; i < e.length; i++) t.children.push(e[i]);
			this.processScene(t)
		},
		processInput: function(e) {
			var t = this.options;
			e = e instanceof Array ? e : [e], this._invokeAll((function(t) {
				t.beforeParse && t.beforeParse(e)
			}));
			for (var i = [], n = 0; n < e.length; n++) e[n] instanceof Scene ? this.processScene(e[n]) : i.push(e[n]);
			i.length > 0 && this.processObjects(i);
			for (n = 0; n < this.skins.length; ++n) this.processSkin(this.skins[n]);
			for (n = 0; n < t.animations.length; ++n) this.processAnimation(t.animations[n], e[0]);
			this._invokeAll((function(t) {
				t.afterParse && t.afterParse(e)
			}))
		},
		_invokeAll: function(e) {
			for (var t = 0, i = this.plugins.length; t < i; t++) e(this.plugins[t])
		}
	}, M.prototype = {
		constructor: M,
		writeNode: function(e, t) {
			if (e.isLight)
				if (e.isDirectionalLight || e.isPointLight || e.isSpotLight) {
					var i = this.writer,
						n = i.json,
						r = i.extensionsUsed,
						a = {};
					e.name && (a.name = e.name), a.color = e.color.toArray(), a.intensity = e.intensity, e.isDirectionalLight ? a.type = "directional" : e.isPointLight ? (a.type = "point", e.distance > 0 && (a.range = e.distance)) : e.isSpotLight && (a.type = "spot", e.distance > 0 && (a.range = e.distance), a.spot = {}, a.spot.innerConeAngle = (e.penumbra - 1) * e.angle * -1, a.spot.outerConeAngle = e.angle), void 0 !== e.decay && 2 !== e.decay && console.warn("THREE.GLTFExporter: Light decay may be lost. glTF is physically-based, and expects light.decay=2."), !e.target || e.target.parent === e && 0 === e.target.position.x && 0 === e.target.position.y && -1 === e.target.position.z || console.warn("THREE.GLTFExporter: Light direction may be lost. For best results, make light.target a child of the light with position 0,0,-1."), r[this.name] || (n.extensions = n.extensions || {}, n.extensions[this.name] = {
						lights: []
					}, r[this.name] = !0);
					var s = n.extensions[this.name].lights;
					s.push(a), t.extensions = t.extensions || {}, t.extensions[this.name] = {
						light: s.length - 1
					}
				} else console.warn("THREE.GLTFExporter: Only directional, point, and spot lights are supported.", e)
		}
	}, T.prototype = {
		constructor: T,
		writeMaterial: function(e, t) {
			if (e.isMeshBasicMaterial) {
				var i = this.writer.extensionsUsed;
				t.extensions = t.extensions || {}, t.extensions[this.name] = {}, i[this.name] = !0, t.pbrMetallicRoughness.metallicFactor = 0, t.pbrMetallicRoughness.roughnessFactor = .9
			}
		}
	}, B.prototype = {
		constructor: B,
		writeMaterial: function(e, t) {
			if (e.isGLTFSpecularGlossinessMaterial) {
				var i = this.writer,
					n = i.extensionsUsed,
					r = {};
				t.pbrMetallicRoughness.baseColorFactor && (r.diffuseFactor = t.pbrMetallicRoughness.baseColorFactor);
				var a = [1, 1, 1];
				if (e.specular.toArray(a, 0), r.specularFactor = a, r.glossinessFactor = e.glossiness, t.pbrMetallicRoughness.baseColorTexture && (r.diffuseTexture = t.pbrMetallicRoughness.baseColorTexture), e.specularMap) {
					var s = {
						index: i.processTexture(e.specularMap)
					};
					i.applyTextureTransform(s, e.specularMap), r.specularGlossinessTexture = s
				}
				t.extensions = t.extensions || {}, t.extensions[this.name] = r, n[this.name] = !0
			}
		}
	}, e.Utils = {
		insertKeyframe: function(e, t) {
			var i, n = e.getValueSize(),
				r = new e.TimeBufferType(e.times.length + 1),
				a = new e.ValueBufferType(e.values.length + n),
				s = e.createInterpolant(new e.ValueBufferType(n));
			if (0 === e.times.length) {
				r[0] = t;
				for (var o = 0; o < n; o++) a[o] = 0;
				i = 0
			} else if (t < e.times[0]) {
				if (Math.abs(e.times[0] - t) < .001) return 0;
				r[0] = t, r.set(e.times, 1), a.set(s.evaluate(t), 0), a.set(e.values, n), i = 0
			} else if (t > e.times[e.times.length - 1]) {
				if (Math.abs(e.times[e.times.length - 1] - t) < .001) return e.times.length - 1;
				r[r.length - 1] = t, r.set(e.times, 0), a.set(e.values, 0), a.set(s.evaluate(t), e.values.length), i = r.length - 1
			} else
				for (o = 0; o < e.times.length; o++) {
					if (Math.abs(e.times[o] - t) < .001) return o;
					if (e.times[o] < t && e.times[o + 1] > t) {
						r.set(e.times.slice(0, o + 1), 0), r[o + 1] = t, r.set(e.times.slice(o + 1), o + 2), a.set(e.values.slice(0, (o + 1) * n), 0), a.set(s.evaluate(t), (o + 1) * n), a.set(e.values.slice((o + 1) * n), (o + 2) * n), i = o + 1;
						break
					}
				}
			return e.times = r, e.values = a, i
		},
		mergeMorphTargetTracks: function(e, t) {
			for (var i = [], n = {}, r = e.tracks, a = 0; a < r.length; ++a) {
				var s = r[a],
					o = PropertyBinding.parseTrackName(s.name),
					l = PropertyBinding.findNode(t, o.nodeName);
				if ("morphTargetInfluences" === o.propertyName && void 0 !== o.propertyIndex) {
					if (s.createInterpolant !== s.InterpolantFactoryMethodDiscrete && s.createInterpolant !== s.InterpolantFactoryMethodLinear) {
						if (s.createInterpolant.isInterpolantFactoryMethodGLTFCubicSpline) throw new Error("THREE.GLTFExporter: Cannot merge tracks with glTF CUBICSPLINE interpolation.");
						console.warn("THREE.GLTFExporter: Morph target interpolation mode not yet supported. Using LINEAR instead."), (s = s.clone()).setInterpolation(2301)
					}
					var c, h = l.morphTargetInfluences.length,
						u = l.morphTargetDictionary[o.propertyIndex];
					if (void 0 === u) throw new Error("THREE.GLTFExporter: Morph target name not found: " + o.propertyIndex);
					if (void 0 !== n[l.uuid]) {
						var d = s.createInterpolant(new s.ValueBufferType(1));
						c = n[l.uuid];
						for (A = 0; A < c.times.length; A++) c.values[A * h + u] = d.evaluate(c.times[A]);
						for (A = 0; A < s.times.length; A++) {
							var p = this.insertKeyframe(c, s.times[A]);
							c.values[p * h + u] = s.values[A]
						}
					} else {
						for (var m = new((c = s.clone()).ValueBufferType)(h * c.times.length), A = 0; A < c.times.length; A++) m[A * h + u] = c.values[A];
						c.name = (o.nodeName || "") + ".morphTargetInfluences", c.values = m, n[l.uuid] = c, i.push(c)
					}
				} else i.push(s)
			}
			return e.tracks = i, e
		}
	}, e
}();
const $correlatedObjects = Symbol("correlatedObjects"),
	$sourceObject = Symbol("sourceObject"),
	$onUpdate = Symbol("onUpdate");
class ThreeDOMElement {
	constructor(e, t, i = null) {
		this[$onUpdate] = e, this[$sourceObject] = t, this[$correlatedObjects] = i
	}
}
var _a$3, _b$2;
const loader = new ImageLoader,
	$threeTextures$1 = Symbol("threeTextures"),
	$uri = Symbol("uri"),
	$bufferViewImages = Symbol("bufferViewImages");
class Image$1 extends ThreeDOMElement {
	constructor(e, t, i) {
		if (super(e, t, i), this[_a$3] = void 0, this[_b$2] = new WeakMap, null != t.uri && (this[$uri] = t.uri), null != t.bufferView)
			for (const e of i) this[$bufferViewImages].set(e, e.image)
	}
	get[$threeTextures$1]() {
		return this[$correlatedObjects]
	}
	get name() {
		return this[$sourceObject].name || ""
	}
	get uri() {
		return this[$uri]
	}
	get type() {
		return null != this.uri ? "external" : "embedded"
	}
	async setURI(e) {
		this[$uri] = e;
		const t = await new Promise((t, i) => {
			loader.load(e, t, void 0, i)
		});
		for (const e of this[$threeTextures$1]) null == t && null != this[$sourceObject].bufferView ? e.image = this[$bufferViewImages].get(e) : e.image = t, e.needsUpdate = !0;
		this[$onUpdate]()
	}
}
_a$3 = $uri, _b$2 = $bufferViewImages;
const isMinFilter = (() => {
		const e = [9728, 9729, 9984, 9985, 9986, 9987];
		return t => e.indexOf(t) > -1
	})(),
	isMagFilter = (() => {
		const e = [9728, 9729];
		return t => e.indexOf(t) > -1
	})(),
	isWrapMode = (() => {
		const e = [33071, 33648, 10497];
		return t => e.indexOf(t) > -1
	})(),
	isValidSamplerValue = (e, t) => {
		switch (e) {
			case "minFilter":
				return isMinFilter(t);
			case "magFilter":
				return isMagFilter(t);
			case "wrapS":
			case "wrapT":
				return isWrapMode(t);
			default:
				throw new Error(`Cannot configure property "${e}" on Sampler`)
		}
	},
	$threeTextures = Symbol("threeTextures"),
	$setProperty = Symbol("setProperty");
class Sampler extends ThreeDOMElement {
	get[$threeTextures]() {
		return this[$correlatedObjects]
	}
	constructor(e, t, i) {
		null == t.minFilter && (t.minFilter = 9987), null == t.magFilter && (t.magFilter = 9729), null == t.wrapS && (t.wrapS = 10497), null == t.wrapT && (t.wrapT = 10497), super(e, t, i)
	}
	get name() {
		return this[$sourceObject].name || ""
	}
	get minFilter() {
		return this[$sourceObject].minFilter
	}
	get magFilter() {
		return this[$sourceObject].magFilter
	}
	get wrapS() {
		return this[$sourceObject].wrapS
	}
	get wrapT() {
		return this[$sourceObject].wrapT
	}
	setMinFilter(e) {
		this[$setProperty]("minFilter", e)
	}
	setMagFilter(e) {
		this[$setProperty]("magFilter", e)
	}
	setWrapS(e) {
		this[$setProperty]("wrapS", e)
	}
	setWrapT(e) {
		this[$setProperty]("wrapT", e)
	} [$setProperty](e, t) {
		const i = this[$sourceObject];
		if (isValidSamplerValue(e, t)) {
			i[e] = t;
			for (const i of this[$threeTextures]) i[e] = t, i.needsUpdate = !0
		}
		this[$onUpdate]()
	}
}
const $source = Symbol("source"),
	$sampler = Symbol("sampler");
class Texture extends ThreeDOMElement {
	constructor(e, t, i, n) {
		super(e, i, n);
		const {
			sampler: r,
			source: a
		} = i, s = null != t.samplers && null != r ? t.samplers[r] : {};
		if (this[$sampler] = new Sampler(e, s, n), null != t.images && null != a) {
			const i = t.images[a];
			null != i && (this[$source] = new Image$1(e, i, n))
		}
	}
	get name() {
		return this[$sourceObject].name || ""
	}
	get sampler() {
		return this[$sampler]
	}
	get source() {
		return this[$source]
	}
}
const $texture = Symbol("texture");
class TextureInfo extends ThreeDOMElement {
	constructor(e, t, i, n) {
		super(e, i, n);
		const {
			index: r
		} = i, a = t.textures[r];
		null != a && (this[$texture] = new Texture(e, t, a, n))
	}
	get texture() {
		return this[$texture]
	}
}
var _a$2, _b$1;
const $threeMaterials = Symbol("threeMaterials"),
	$baseColorTexture = Symbol("baseColorTexture"),
	$metallicRoughnessTexture = Symbol("metallicRoughnessTexture");
class PBRMetallicRoughness extends ThreeDOMElement {
	constructor(e, t, i, n) {
		super(e, i, n), this[_a$2] = null, this[_b$1] = null, null == i.baseColorFactor && (i.baseColorFactor = [1, 1, 1, 1]), null == i.roughnessFactor && (i.roughnessFactor = 0), null == i.metallicFactor && (i.metallicFactor = 0);
		const {
			baseColorTexture: r,
			metallicRoughnessTexture: a
		} = i, s = new Set, o = new Set;
		for (const e of n) null != r && null != e.map && s.add(e.map), null != a && null != e.metalnessMap && o.add(e.metalnessMap);
		s.size > 0 && (this[$baseColorTexture] = new TextureInfo(e, t, r, s)), o.size > 0 && (this[$metallicRoughnessTexture] = new TextureInfo(e, t, a, o))
	}
	get[(_a$2 = $baseColorTexture, _b$1 = $metallicRoughnessTexture, $threeMaterials)]() {
		return this[$correlatedObjects]
	}
	get baseColorFactor() {
		return this[$sourceObject].baseColorFactor
	}
	get metallicFactor() {
		return this[$sourceObject].metallicFactor
	}
	get roughnessFactor() {
		return this[$sourceObject].roughnessFactor
	}
	get baseColorTexture() {
		return this[$baseColorTexture]
	}
	get metallicRoughnessTexture() {
		return this[$metallicRoughnessTexture]
	}
	setBaseColorFactor(e) {
		for (const t of this[$threeMaterials]) t.color.fromArray(e), t.opacity = e[3];
		this[$sourceObject].baseColorFactor = e, this[$onUpdate]()
	}
	setMetallicFactor(e) {
		for (const t of this[$threeMaterials]) t.metalness = e;
		this[$sourceObject].metallicFactor = e, this[$onUpdate]()
	}
	setRoughnessFactor(e) {
		for (const t of this[$threeMaterials]) t.roughness = e;
		this[$sourceObject].roughnessFactor = e, this[$onUpdate]()
	}
}
var _a$1, _b, _c;
const $pbrMetallicRoughness = Symbol("pbrMetallicRoughness"),
	$normalTexture = Symbol("normalTexture"),
	$occlusionTexture = Symbol("occlusionTexture"),
	$emissiveTexture = Symbol("emissiveTexture");
class Material extends ThreeDOMElement {
	constructor(e, t, i, n) {
		if (super(e, i, n), this[_a$1] = null, this[_b] = null, this[_c] = null, null == n) return;
		null == i.pbrMetallicRoughness && (i.pbrMetallicRoughness = {}), this[$pbrMetallicRoughness] = new PBRMetallicRoughness(e, t, i.pbrMetallicRoughness, n);
		const {
			normalTexture: r,
			occlusionTexture: a,
			emissiveTexture: s
		} = i, o = new Set, l = new Set, c = new Set;
		for (const e of n) {
			const {
				normalMap: t,
				aoMap: i,
				emissiveMap: n
			} = e;
			null != r && null != t && o.add(t), null != a && null != i && l.add(i), null != s && null != n && c.add(n)
		}
		o.size > 0 && (this[$normalTexture] = new TextureInfo(e, t, r, o)), l.size > 0 && (this[$occlusionTexture] = new TextureInfo(e, t, a, l)), c.size > 0 && (this[$emissiveTexture] = new TextureInfo(e, t, s, c))
	}
	get name() {
		return this[$sourceObject].name || ""
	}
	get pbrMetallicRoughness() {
		return this[$pbrMetallicRoughness]
	}
	get normalTexture() {
		return this[$normalTexture]
	}
	get occlusionTexture() {
		return this[$occlusionTexture]
	}
	get emissiveTexture() {
		return this[$emissiveTexture]
	}
	get emissiveFactor() {
		return this[$sourceObject].emissiveFactor
	}
	setEmissiveFactor(e) {
		for (const t of this[$correlatedObjects]) t.emissive.fromArray(e);
		this[$sourceObject].emissiveFactor = e, this[$onUpdate]()
	}
}
var _a;
_a$1 = $normalTexture, _b = $occlusionTexture, _c = $emissiveTexture;
const $materials = Symbol("materials");
class Model {
	constructor(e, t = (() => {})) {
		this[_a] = [];
		const {
			gltf: i,
			gltfElementMap: n
		} = e;
		i.materials.forEach(e => {
			this[$materials].push(new Material(t, i, e, n.get(e)))
		})
	}
	get materials() {
		return this[$materials]
	}
}
_a = $materials;
var __decorate$1 = function(e, t, i, n) {
	for (var r, a = arguments.length, s = a < 3 ? t : null === n ? n = Object.getOwnPropertyDescriptor(t, i) : n, o = e.length - 1; o >= 0; o--)(r = e[o]) && (s = (a < 3 ? r(s) : a > 3 ? r(t, i, s) : r(t, i)) || s);
	return a > 3 && s && Object.defineProperty(t, i, s), s
};
const $currentGLTF = Symbol("currentGLTF"),
	$model = Symbol("model"),
	$variants = Symbol("variants"),
	SceneGraphMixin = e => {
		var t, i, n;
		class r extends e {
			constructor() {
				super(...arguments), this[t] = void 0, this[i] = null, this[n] = [], this.variantName = void 0, this.orientation = "0 0 0", this.scale = "1 1 1"
			}
			get model() {
				return this[$model]
			}
			get availableVariants() {
				return this[$variants]
			}
			updated(e) {
				if (super.updated(e), e.has("variantName")) {
					const e = this[$variants],
						t = this[$currentGLTF],
						{
							variantName: i
						} = this,
						n = e.findIndex(e => e === i);
					if (null == t || n < 0) return;
					const r = () => {
							this[$needsRender]()
						},
						a = t.correlatedSceneGraph.loadVariant(n, r),
						{
							gltf: s,
							gltfElementMap: o
						} = t.correlatedSceneGraph;
					for (const e of a) {
						const t = s.materials[e];
						this[$model].materials[e] = new Material(r, s, t, o.get(t))
					}
				}
				if (e.has("orientation") || e.has("scale")) {
					const {
						modelContainer: e
					} = this[$scene], t = parseExpressions(this.orientation)[0].terms, i = normalizeUnit(t[0]).number, n = normalizeUnit(t[1]).number, r = normalizeUnit(t[2]).number;
					e.quaternion.setFromEuler(new Euler(n, r, i, "YXZ"));
					const a = parseExpressions(this.scale)[0].terms;
					e.scale.set(a[0].number, a[1].number, a[2].number), this[$scene].updateBoundingBox(), this[$scene].updateShadow(), this[$renderer].arRenderer.onUpdateScene(), this[$needsRender]()
				}
			} [(t = $model, i = $currentGLTF, n = $variants, $onModelLoad)]() {
				super[$onModelLoad](), this[$variants] = [];
				const {
					currentGLTF: e
				} = this[$scene];
				if (null != e) {
					const {
						correlatedSceneGraph: t
					} = e;
					null != t && e !== this[$currentGLTF] && (this[$model] = new Model(t, () => {
						this[$needsRender]()
					}));
					const {
						gltfExtensions: i
					} = e.userData;
					if (null != i) {
						const e = i.KHR_materials_variants;
						null != e && (this[$variants] = e.variants.map(e => e.name), this.requestUpdate("variantName"))
					}
				}
				this[$currentGLTF] = e, this.dispatchEvent(new CustomEvent("scene-graph-ready"))
			}
			async exportScene(e) {
				const t = this[$scene];
				return new Promise(async i => {
					const n = {
						binary: !0,
						onlyVisible: !0,
						maxTextureSize: 1 / 0,
						forcePowerOfTwoTextures: !1,
						includeCustomExtensions: !1,
						embedImages: !0
					};
					Object.assign(n, e), n.animations = t.animations, n.truncateDrawRange = !0;
					const r = t.shadow;
					let a = !1;
					null != r && (a = r.visible, r.visible = !1);
					(new GLTFExporter).parse(t.modelContainer, e => i(new Blob([n.binary ? e : JSON.stringify(e)], {
						type: n.binary ? "application/octet-stream" : "application/json"
					})), n), null != r && (r.visible = a)
				})
			}
		}
		return __decorate$1([property({
			type: String,
			attribute: "variant-name"
		})], r.prototype, "variantName", void 0), __decorate$1([property({
			type: String,
			attribute: "orientation"
		})], r.prototype, "orientation", void 0), __decorate$1([property({
			type: String,
			attribute: "scale"
		})], r.prototype, "scale", void 0), r
	};
var __decorate = function(e, t, i, n) {
	for (var r, a = arguments.length, s = a < 3 ? t : null === n ? n = Object.getOwnPropertyDescriptor(t, i) : n, o = e.length - 1; o >= 0; o--)(r = e[o]) && (s = (a < 3 ? r(s) : a > 3 ? r(t, i, s) : r(t, i)) || s);
	return a > 3 && s && Object.defineProperty(t, i, s), s
};
const DEFAULT_ROTATION_SPEED = Math.PI / 32,
	AUTO_ROTATE_DELAY_DEFAULT = 3e3,
	rotationRateIntrinsics = {
		basis: [degreesToRadians(numberNode(DEFAULT_ROTATION_SPEED, "rad"))],
		keywords: {
			auto: [null]
		}
	},
	$autoRotateStartTime = Symbol("autoRotateStartTime"),
	$radiansPerSecond = Symbol("radiansPerSecond"),
	$syncRotationRate = Symbol("syncRotationRate"),
	$onCameraChange = Symbol("onCameraChange"),
	StagingMixin = e => {
		var t, i, n;
		class r extends e {
			constructor() {
				super(...arguments), this.autoRotate = !1, this.autoRotateDelay = 3e3, this.rotationPerSecond = "auto", this[t] = performance.now(), this[i] = 0, this[n] = e => {
					this.autoRotate && "user-interaction" === e.detail.source && (this[$autoRotateStartTime] = performance.now())
				}
			}
			connectedCallback() {
				super.connectedCallback(), this.addEventListener("camera-change", this[$onCameraChange]), this[$autoRotateStartTime] = performance.now()
			}
			disconnectedCallback() {
				super.disconnectedCallback(), this.removeEventListener("camera-change", this[$onCameraChange]), this[$autoRotateStartTime] = performance.now()
			}
			updated(e) {
				super.updated(e), e.has("autoRotate") && (this[$autoRotateStartTime] = performance.now())
			} [(t = $autoRotateStartTime, i = $radiansPerSecond, $syncRotationRate)](e) {
				this[$radiansPerSecond] = e[0]
			} [$tick](e, t) {
				if (super[$tick](e, t), !this.autoRotate || !this[$hasTransitioned]() || this[$renderer].isPresenting) return;
				const i = Math.min(t, e - this[$autoRotateStartTime] - this.autoRotateDelay);
				i > 0 && (this[$scene].yaw = this.turntableRotation + this[$radiansPerSecond] * i * .001)
			}
			get turntableRotation() {
				return this[$scene].yaw
			}
			resetTurntableRotation(e = 0) {
				this[$scene].yaw = e
			}
		}
		return n = $onCameraChange, __decorate([property({
			type: Boolean,
			attribute: "auto-rotate"
		})], r.prototype, "autoRotate", void 0), __decorate([property({
			type: Number,
			attribute: "auto-rotate-delay"
		})], r.prototype, "autoRotateDelay", void 0), __decorate([style({
			intrinsics: rotationRateIntrinsics,
			updateHandler: $syncRotationRate
		}), property({
			type: String,
			attribute: "rotation-per-second"
		})], r.prototype, "rotationPerSecond", void 0), r
	},
	FocusVisiblePolyfillMixin = e => {
		var t;
		const i = Symbol("endPolyfillCoordination");
		return t = i, class extends e {
			constructor() {
				super(...arguments), this[t] = null
			}
			connectedCallback() {
				super.connectedCallback && super.connectedCallback(), null == this[i] && (this[i] = (e => {
					if (null == e.shadowRoot || e.hasAttribute("data-js-focus-visible")) return () => {};
					if (!self.applyFocusVisiblePolyfill) {
						const t = () => {
							self.applyFocusVisiblePolyfill(e.shadowRoot)
						};
						return self.addEventListener("focus-visible-polyfill-ready", t, {
							once: !0
						}), () => {
							self.removeEventListener("focus-visible-polyfill-ready", t)
						}
					}
					return self.applyFocusVisiblePolyfill(e.shadowRoot), () => {}
				})(this))
			}
			disconnectedCallback() {
				super.disconnectedCallback && super.disconnectedCallback(), null != this[i] && (this[i](), this[i] = null)
			}
		}
	},
	ModelViewerElement = AnnotationMixin(SceneGraphMixin(StagingMixin(EnvironmentMixin(ControlsMixin(ARMixin(LoadingMixin(AnimationMixin(FocusVisiblePolyfillMixin(ModelViewerElementBase)))))))));
customElements.define("model-viewer", ModelViewerElement);
const body = document.querySelector("body");
body.classList.add("--lock"), window.addEventListener("DOMContentLoaded", () => {
	const e = document.querySelector(".loader");
	setTimeout(() => {
		e.style.opacity = "0", e.style.transform = "scale(2) ", setTimeout(() => {
			e.style.display = "none"
		}, 500)
	}, 1e3);
	const t = document.querySelectorAll(".--show-item");
	if (t.length > 0)
		for (let e = 0; e < t.length; e++) {
			const i = t[e];
			i.classList.contains("--show") || i.classList.add("--show")
		}
}), setTimeout(() => {
	body.classList.remove("--lock");
	const e = document.querySelector(".about");
	e.classList.contains("--bgAnimate") || e.classList.add("--bgAnimate")
}, 9500);
const header = document.querySelector(".header"),
	burgerMenu = document.querySelector(".burger-menu"),
	navbar = document.querySelector(".menu__navbar"),
	listItems = document.querySelectorAll(".menu__list-item"),
	listItem1 = document.querySelector(".menu__list-item-1"),
	listItem2 = document.querySelector(".menu__list-item-2"),
	listItem3 = document.querySelector(".menu__list-item-3"),
	listItem4 = document.querySelector(".menu__list-item-4"),
	langLink = document.querySelector(".menu__link-4"),
	langMenu = document.querySelector(".menu__lang-list");
burgerMenu.addEventListener("click", (function(e) {
	body.classList.toggle("--lock--mobile"), header.classList.toggle("--active"), burgerMenu.classList.toggle("--active"), navbar.classList.toggle("--active");
	for (let e = 0; e < listItems.length; e++) {
		listItems[e].classList.toggle("--active")
	}
	langMenu.classList.remove("--active")
}));
for (let e = 0; e < listItems.length; e++) {
	const t = listItems[e];
	t !== listItems[3] && t.addEventListener("click", (function() {
		body.classList.remove("--lock--modbile"), header.classList.remove("--active"), burgerMenu.classList.remove("--active"), navbar.classList.remove("--active"), listItem1.classList.remove("--active"), listItem2.classList.remove("--active"), listItem3.classList.remove("--active"), listItem4.classList.remove("--active"), langMenu.classList.remove("--active")
	}))
}
langLink.addEventListener("click", (function(e) {
	langMenu.classList.toggle("--active"), langMenu.addEventListener("click", (function(e) {
		e.target.closest("menu__lang-list") || langMenu.classList.remove("--active")
	}))
}));
const swiper = new Swiper(".main-content__slider", {
		loop: !0,
		parallax: !0,
		pagination: {
			el: ".swiper-pagination",
			clickable: !0,
			dynamicBullets: !0
		},
		navigation: {
			nextEl: ".swiper-button-next",
			prevEl: ".swiper-button-prev"
		},
		autoplay: {
			delay: 3e3,
			stopOnLastSlide: !1,
			disableOnInteraction: !1
		},
		speed: 1500,
		grabCursor: !0,
		slideToClickedSlide: !0,
		keyboard: {
			enabled: !0,
			onlyInViewport: !0,
			pageUpDown: !0
		},
		watchOverflow: !0,
		effect: "fade",
		fadeEffect: {
			crossFade: !0
		},
		preloadImages: !1,
		watchSlidesProgress: !0,
		watchSlidesVisibility: !0,
		a11y: {
			prevSlideMessage: "Previous slide",
			nextSlideMessage: "Next slide",
			firstSlideMessage: "This is the first slide",
			lastSlideMessage: "This is the first slide",
			notificationClass: "swiper-notification",
			paginationBulletMessage: "Go to slide {{index}}"
		}
	}),
	popupLinks = document.querySelectorAll(".popup-link"),
	lockPadding = document.querySelectorAll(".lock-padding");
let unlock = !0;
const timeout = 200;
if (popupLinks.length > 0)
	for (let e = 0; e < popupLinks.length; e++) {
		const t = popupLinks[e];
		t.addEventListener("click", (function(e) {
			const i = t.getAttribute("href").replace("#", "");
			popupOpen(document.getElementById(i)), e.preventDefault()
		}))
	}
const popupCloseIcon = document.querySelectorAll(".close-popup");
if (popupCloseIcon.length > 0)
	for (let e = 0; e < popupCloseIcon.length; e++) {
		const t = popupCloseIcon[e];
		t.addEventListener("click", (function(e) {
			popupClose(t.closest(".popup")), e.preventDefault()
		}))
	}

function popupOpen(e) {
	if (e && unlock) {
		const t = document.querySelector(".popup.--open");
		t ? popupClose(t, !1) : bodyLock(), e.classList.add("--open"), e.addEventListener("click", (function(e) {
			e.target.closest(".popup__content") || popupClose(e.target.closest(".popup"))
		}))
	}
}

function popupClose(e, t = !0) {
	unlock && (e.classList.remove("--open"), t && bodyUnlock())
}

function bodyLock() {
	const e = window.innerWidth - document.querySelector(".wrapper").offsetWidth + "px";
	if (lockPadding.length > 0)
		for (let t = 0; t < lockPadding.length; t++) {
			lockPadding[t].style.paddingRight = e
		}
	body.style.paddingRight = e, body.classList.add("--lock"), unlock = !1, setTimeout((function() {
		unlock = !0
	}), 200)
}

function bodyUnlock() {
	setTimeout(() => {
		if (lockPadding.length > 0)
			for (let e = 0; e < lockPadding.length; e++) {
				lockPadding[e].style.paddingRight = "0px"
			}
		body.style.paddingRight = "0px", body.classList.remove("--lock")
	}, 200), unlock = !1, setTimeout(() => {
		unlock = !0
	}, 200)
}
document.addEventListener("keydown", (function(e) {
	if (27 === e.which) {
		const e = document.querySelectorAll(".popup.--open");
		if (e.length > 0)
			for (let t = 0; t < e.length; t++) {
				popupClose(e[t])
			}
	}
})), document.querySelector(".popup").classList.contains("--open") && body.classList.add("--lock");
const aboutLink = document.querySelector(".about__link");
aboutLink.addEventListener("mousedown", (function() {
	aboutLink.style.transform = "scale(0.9)", aboutLink.addEventListener("mousedown", (function() {
		aboutLink.style.transform = "scale(1)"
	}))
}));
const contactButton = document.querySelector(".form__button");
contactButton.addEventListener("click", (function() {
	clientName(), clientNumber()
}));
const form = document.querySelector("#form");

function clientName() {
	let e, t, i = document.getElementById("firstName");
	t = document.getElementById("contact__error-1"), e = document.getElementById("firstName").value;
	try {
		if ("" == e) throw i.classList.add("contact--error"), i.classList.remove("contact--valid"), "!";
		e.length >= 3 && (i.classList.remove("contact--error"), i.classList.add("contact--valid"), t.innerHTML = "")
	} catch (e) {
		t.classList.contains("--ru") ? t.innerHTML = "Введите свое имя" + e : t.classList.contains("--en") ? t.innerHTML = "Enter your name" + e : t.innerHTML = "Ismingizni kiriting" + e
	}
}

function clientNumber() {
	let e, t, i = document.getElementById("phoneNumber");
	t = document.getElementById("contact__error-2"), e = document.getElementById("phoneNumber").value;
	try {
		if ("" == e) throw i.classList.add("contact--error"), i.classList.remove("contact--valid"), "!";
		13 == e.length ? (i.classList.remove("contact--error"), i.classList.add("contact--valid"), t.innerHTML = "") : (t.classList.contains("--ru") ? t.innerHTML = "Введите свой номер в виде +998 12 345 67 89!" : t.classList.contains("--en") ? t.innerHTML = "Enter your number as  +998 12 345 67 89!" : t.innerHTML = "Raqamingizni  +998 12 345 67 89 ko'rinishda kiriting!", i.classList.add("contact--error"), i.classList.remove("contact--valid"))
	} catch (e) {
		t.classList.contains("--ru") ? t.innerHTML = "Введите свой номер" + e : t.classList.contains("--en") ? t.innerHTML = "Enter your number" + e : t.innerHTML = "Raqamingizni  kiriting" + e
	}
}

function isValueName(e) {
	let t = String.fromCharCode(e.which);
	/[0-9-+-@-!-_-#-$-%-^-&-*]/.test(t) && !/[']/.test(t) && e.preventDefault()
}

function isValueNum(e) {
	let t = String.fromCharCode(e.which);
	/[0-9]/.test(t) || e.preventDefault()
}
form.addEventListener("submit", e => {
	e.preventDefault();
	const t = `https://api.telegram.org/bot1811116830:AAEMSCDonDOQZhEct6VtAiJ4vUBOFStsRxs/sendMessage?chat_id=-346912975&text=${`Eropen.uz <b>Online ariza</b> 📬 %0A----------------------------------%0A <b>👤Ismi: </b> <i>${document.getElementById("firstName").value}</i>%0A  <b>📞Tel.raqami: /[+]/ </b> <i>${document.getElementById("phoneNumber").value}</i>%0A <b>✉️Text: </b> <i>${document.getElementById("message").value}</i>`}&parse_mode=html`;
	let i = new XMLHttpRequest;
	i.open("GET", t, !0), i.send();
	const n = document.querySelector("#contact__popup"),
		r = document.querySelector("#contact__popup-error");
	document.getElementById("firstName").value.length >= 3 && 13 == document.getElementById("phoneNumber").value.length && document.getElementById("message").value.length > 0 ? (n.classList.add("--open"), body.classList.add("--lock"), document.querySelector(".contact__user-name").innerHTML = document.getElementById("firstName").value, document.querySelector(".contact__user-number").innerHTML = document.getElementById("phoneNumber").value, r.classList.remove("--open"), form.reset(), document.getElementById("firstName").classList.remove("contact--valid"), document.getElementById("phoneNumber").classList.remove("contact--valid"), setTimeout(() => {
		body.classList.remove("--lock"), n.classList.remove("--open")
	}, 8e3)) : (n.classList.remove("--open"), r.classList.add("--open"))
}), contactButton.addEventListener("mousedown", (function() {
	contactButton.style.transform = "scale(0.9)", contactButton.addEventListener("mouseup", (function() {
		contactButton.style.transform = "scale(1)"
	}))
})), document.getElementById("phoneNumber").addEventListener("keyup", (function() {
	clientNumber()
})), document.getElementById("firstName").addEventListener("keyup", (function() {
	clientName()
}));
const animItems = document.querySelectorAll(".--anim-items");
if (animItems.length > 0) {
	function animOnScroll() {
		for (let t = 0; t < animItems.length; t++) {
			const i = animItems[t],
				n = i.offsetHeight,
				r = e(i).top,
				a = 6;
			let s = window.innerHeight - n / a;
			n > window.innerHeight && (s = window.innerHeight - window.innerHeight / a), pageYOffset > r - s && pageYOffset < r + n ? i.classList.add("--animate") : i.classList.remove("--animate")
		}

		function e(e) {
			const t = e.getBoundingClientRect(),
				i = window.pageXOffset || document.documentElement.scrollLeft;
			return scrollTop = window.pageYOffset || document.documentElement.scrollTop, {
				top: t.top + scrollTop,
				left: t.left + i
			}
		}
	}
	window.addEventListener("scroll", animOnScroll), setTimeout(() => {
		animOnScroll()
	}, 400)
}
const collapseButton = document.querySelector(".collapse-button");
collapseButton.addEventListener("click", (function() {
	collapseButton.nextElementSibling.hidden = !0
}));
const spollersArray = document.querySelectorAll("[data-spollers]");
if (spollersArray.length > 0) {
	const e = Array.from(spollersArray).filter((function(e, t, i) {
		return !e.dataset.spollers.split(",")[0]
	}));

	function initSpollers(e, t = !1) {
		e.forEach(e => {
			t ? (e.classList.remove("--init"), initSpollersBody(e, !1), e.removeEventListener("click", setSpollerAction)) : (e.classList.add("--init"), initSpollersBody(e), e.addEventListener("click", setSpollerAction))
		})
	}

	function initSpollersBody(e, t = !0) {
		const i = e.querySelectorAll("[data-spoller]");
		i.length > 0 && i.forEach(e => {
			t ? (e.removeAttribute("tabindex"), e.classList.contains("--active") || (e.nextElementSibling.hidden = !0)) : (e.setAttribute("tabindex", "-1"), e.nextElementSibling.hidden = !1)
		})
	}

	function setSpollerAction(e) {
		const t = e.target;
		if (t.hasAttribute("data-spoller") || t.closest("[data-spoller]")) {
			const i = t.hasAttribute("data-spoller") ? t : t.closest("[data-spoller]"),
				n = i.closest("[data-spollers]"),
				r = !!n.hasAttribute("data-one-spoller");
			n.querySelectorAll(".--slide").length || (r && !i.classList.contains("--active") && hideSpollersBody(n), i.classList.toggle("--active"), slideToggle(i.nextElementSibling, 500)), e.preventDefault()
		}
	}

	function hideSpollersBody(e) {
		const t = e.querySelector("[data-spoller].--active");
		t && (t.classList.remove("--active"), slideUp(t.nextElementSibling, 500))
	}
	e.length > 0 && initSpollers(e)
}
let slideUp = (e, t = 500) => {
		e.classList.contains("--slide") || (e.classList.add("--slide"), e.style.transitionProperty = "height, margin, padding", e.style.transitionDuration = t + "ms", e.style.height = e.offsetHeight + "px", e.offsetHeight, e.style.overflow = "hidden", e.style.height = 0, e.style.paddingTop = 0, e.style.paddingBottom = 0, e.style.marginTop = 0, e.style.marginBottom = 0, window.setTimeout(() => {
			e.hidden = !0, e.style.removeProperty("height"), e.style.removeProperty("padding-top"), e.style.removeProperty("padding-bottom"), e.style.removeProperty("margin-top"), e.style.removeProperty("margin-bottom"), e.style.removeProperty("overflow"), e.style.removeProperty("transition-duration"), e.style.removeProperty("transition-property"), e.classList.remove("--slide")
		}, t))
	},
	slideDown = (e, t = 500) => {
		if (!e.classList.contains("--slide")) {
			e.classList.add("--slide"), e.hidden && (e.hidden = !1);
			let i = e.offsetHeight;
			e.style.overflow = "hidden", e.style.height = 0, e.style.paddingTop = 0, e.style.paddingBottom = 0, e.style.marginTop = 0, e.style.marginBottom = 0, e.offsetHeight, e.style.transitionProperty = "height, margin, padding", e.style.transitionDuration = t + "ms", e.style.height = i + "px", e.style.removeProperty("padding-top"), e.style.removeProperty("padding-bottom"), e.style.removeProperty("margin-top"), e.style.removeProperty("margin-bottom"), window.setTimeout(() => {
				e.style.removeProperty("height"), e.style.removeProperty("overflow"), e.style.removeProperty("transition-duration"), e.style.removeProperty("transition-property"), e.classList.remove("--slide")
			}, t)
		}
	},
	slideToggle = (e, t = 500) => e.hidden ? slideDown(e, t) : slideUp(e, t);