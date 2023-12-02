@tool
class_name OctasphereMesh extends ArrayMesh

@export_range(0, 8) var subdivisions : int = 2:
	set(value):
		if not value == subdivisions:
			subdivisions = value
			update_mesh()

@export var radius := 1.0 :
	set(value):
		if not radius ==  value:
			radius = value
			update_mesh()

@export var width := 0.0 :
	set(value):
		if not width ==  value:
			width = value
			update_mesh()
			
@export var height := 0.0 :
	set(value):
		if not height ==  value:
			height = value
			update_mesh()
			
@export var depth := 0.0 :
	set(value):
		if not depth ==  value:
			depth = value
			update_mesh()
			
func _init():
	if self.get_surface_count() == 0:
		update_mesh()

func update_mesh():
	var ans = octasphere(subdivisions, radius, 0, 0, 0)
	var normals: PackedVector3Array = ans[0]
	ans = octasphere(subdivisions, radius, width, height, depth)
	var vertices: PackedVector3Array = ans[0]
	var triangles: PackedInt32Array = ans[1]
	var tex_uv: PackedVector2Array = ans[2]
	var colors: PackedColorArray = ans[3]

	## Initialize the ArrayMesh.
	var arrays = []
	arrays.resize(ArrayMesh.ARRAY_MAX)
	arrays[ArrayMesh.ARRAY_VERTEX] = vertices
	arrays[ArrayMesh.ARRAY_TEX_UV] = tex_uv
	arrays[ArrayMesh.ARRAY_INDEX] = triangles
	arrays[ArrayMesh.ARRAY_NORMAL] = normals
	arrays[ArrayMesh.ARRAY_COLOR] = colors
	## Create the Mesh.
	clear_surfaces()
	self.add_surface_from_arrays(Mesh.PRIMITIVE_TRIANGLES, arrays)
	emit_changed()

func octasphere(ndivisions: int, radius: float, width=0, height=0, depth=0) -> Array:
	var r2 = 2 * radius
	width = max(width, r2)
	height = max(height, r2)
	depth = max(depth, r2)
	var n = 2**ndivisions + 1
	var num_verts: int = n * (n + 1) / 2
	var verts: PackedVector3Array = []
	verts.resize(num_verts)
	var j = 0
	for i in range(n):
		var theta = PI * 0.5 * i / (n - 1)
		var point_a = Vector3(0, sin(theta), cos(theta))
		var point_b = Vector3(cos(theta), sin(theta), 0)
		var num_segments = n - 1 - i
		j = compute_geodesic(verts, j, point_a, point_b, num_segments)
	assert(verts.size() == num_verts)
	for i in range(num_verts):
		verts[i] *= radius

	var num_faces = (n - 2) * (n - 1) + n - 1
	var faces : PackedVector3Array = []
	var f = 0
	var j0 = 0
	for col_index in range(n-1):
		var col_height = n - 1 - col_index
		var j1 = j0 + 1
		var j2 = j0 + col_height + 1
		var j3 = j0 + col_height + 2
		for row in range(col_height - 1):
			faces.append(Vector3(j0 + row, j1 + row, j2 + row))
			faces.append(Vector3(j2 + row, j1 + row, j3 + row))
			f += 2
		var row = col_height - 1
		faces.append(Vector3(j0 + row, j1 + row, j2 + row))
		f += 1
		j0 = j2

	var euler_angles = [
		Vector3(0, 0, 0), Vector3(0, 1, 0), Vector3(0, 2, 0), Vector3(0, 3, 0),
		Vector3(1, 0, 0), Vector3(1, 0, 1), Vector3(1, 0, 2), Vector3(1, 0, 3),
	]
	var vertex_colors = [
		Color.RED, Color.GREEN, Color.BLUE, Color.YELLOW,
		Color.MAGENTA, Color.CYAN, Color.DARK_KHAKI, Color.BROWN, 
	]
	var quats = []
	for e in euler_angles:
		quats.append(Quaternion.from_euler(e * PI * 0.5))

	var offset = 0
	var combined_verts: PackedVector3Array = []
	var combined_faces: PackedInt32Array = []
	var combined_colors: PackedColorArray = []
	for quat in quats:
		var rotated_verts: PackedVector3Array = []
		for v in verts:
			rotated_verts.append(quat * v)
			combined_colors.append(vertex_colors[offset / num_verts])
		var rotated_faces: PackedInt32Array = []
		for fa in faces:
			rotated_faces.append(fa[2] + offset)
			rotated_faces.append(fa[1] + offset)
			rotated_faces.append(fa[0] + offset)

		combined_verts.append_array(rotated_verts)
		combined_faces.append_array(rotated_faces)
		offset += verts.size()

	var tx: float = (width - r2) / 2.0
	var ty: float = (height - r2) / 2.0
	var tz: float = (depth - r2) / 2.0
	var translation: Vector3 = Vector3(tx, ty, tz)

	if not translation.is_equal_approx(Vector3.ZERO):
		var translation_matrix: Array = [
			Vector3(1, 1, 1), Vector3(1, 1, -1), Vector3(-1, 1, -1), Vector3(-1, 1, 1),
			Vector3(1, -1, 1), Vector3(-1, -1, 1), Vector3(-1, -1, -1), Vector3(1, -1, -1)
		]

		for i in range(translation_matrix.size()):
			translation_matrix[i] *= translation

		for i in range(0, combined_verts.size(), num_verts):
			for ii in range(num_verts):
				if i + ii < combined_verts.size():
					combined_verts[i + ii] += translation_matrix[i / num_verts]

	var connectors = add_connectors(ndivisions, radius, width, height, depth)
	if radius == 0:
		assert(len(connectors) / 2 == 6)
		combined_faces = connectors
	else:
		combined_faces.append_array(connectors)

	var combined_uv: PackedVector2Array = []
	for i in range(combined_verts.size()):
		var octant = i / num_verts
		var relative_index = i % num_verts
		var uv := Vector2(0,0)
		var xyz = combined_verts[i]
		var x = xyz[0]
		var y = xyz[1]
		var z = xyz[2]
		var phi = -atan2(z, x)
		var theta = acos(y)
		uv[0] = 0.5 * (phi / PI + 1.0)
		uv[1] = theta / PI
	
		# Special case for the north pole
		if octant < 4 and relative_index == num_verts - 1:
			uv[0] = fmod(0.375 + 0.25 * octant, 1.0)
			uv[1] = 0
	
		# Special case for the south pole
		if octant >= 4 and relative_index == 0:
			uv[0] = 0.375 - 0.25 * (octant - 4)
			if uv[0] + uv[0] < 0:
				uv[0] += 1.0
			uv[1] = 1.0
	
		# Adjust the prime meridian for proper wrapping
		if (octant == 2 || octant == 6) and uv[0] < 0.5:
			uv[0] += 1.0
		
		combined_uv.append(uv)
		
	return [combined_verts, combined_faces, combined_uv, combined_colors]


func compute_geodesic(dst, index, point_a, point_b, num_segments) -> int:
	var angle_between_endpoints = acos(point_a.dot(point_b))
	var rotation_axis = (point_a.cross(point_b)).normalized()
	dst[index] = point_a
	index += 1
	if num_segments == 0:
		return index
	var dtheta = angle_between_endpoints / num_segments
	for point_index in range(1, num_segments):
		var theta = point_index * dtheta
		var q = Quaternion(rotation_axis, theta)
		dst[index] = q * point_a
		index += 1
	dst[index] = point_b
	return index + 1

func add_connectors(ndivisions: int, radius: float, width: float, height: float, depth: float) -> Array:
	var r2 = 2 * radius
	width = max(width, r2)
	height = max(height, r2)
	depth = max(depth, r2)
	var n = 2**ndivisions + 1
	var num_verts: int = n * (n + 1) / 2
	var tx = (width - r2) / 2
	var ty = (height - r2) / 2
	var tz = (depth - r2) / 2

	var boundaries: Array[PackedInt32Array] = get_boundary_indices(ndivisions)
	assert(boundaries.size() == 3)
	var connectors: PackedInt32Array = []

	if radius > 0:
		# Top half
		for patch in range(4):
			if patch % 2 == 0 and tz == 0:
				continue
			if patch % 2 == 1 and tx == 0:
				continue
			var next_patch = (patch + 1) % 4
			var boundary_a = boundaries[1].duplicate()
			for i in range(boundary_a.size()):
				boundary_a[i] += num_verts * patch
			var boundary_b = boundaries[0].duplicate()
			for i in range(boundary_b.size()):
				boundary_b[i] += num_verts * next_patch
			for i in range(n-1):
				var a = boundary_a[i]
				var b = boundary_b[i]
				var c = boundary_a[i+1]
				var d = boundary_b[i+1]
				connectors.append_array([a, b, d, d, c, a])
		
		# Bottom half
		for patch in range(4, 8):
			if patch % 2 == 0 and tx == 0:
				continue
			if patch % 2 == 1 and tz == 0:
				continue
			var next_patch = 4 + (patch + 1) % 4
			var boundary_a = boundaries[0].duplicate()
			for i in range(boundary_a.size()):
				boundary_a[i] += num_verts * patch
			var boundary_b = boundaries[2].duplicate()
			for i in range(boundary_b.size()):
				boundary_b[i] += num_verts * next_patch
			for i in range(n-1):
				var a = boundary_a[i]
				var b = boundary_b[i]
				var c = boundary_a[i+1]
				var d = boundary_b[i+1]
				connectors.append_array([d, b, a, a, c, d])
		
		# Connect top patch to bottom patch
		if ty > 0:
			for patch in range(4):
				var next_patch = 4 + (4 - patch) % 4
				var boundary_a = boundaries[2].duplicate()
				for i in range(boundary_a.size()):
					boundary_a[i] += num_verts * patch
				var boundary_b = boundaries[1].duplicate()
				for i in range(boundary_b.size()):
					boundary_b[i] += num_verts * next_patch
				for i in range(n-1):
					var a = boundary_a[i]
					var b = boundary_b[n-1-i]
					var c = boundary_a[i+1]
					var d = boundary_b[n-1-i-1]
					connectors.append_array([a, b, d, d, c, a])

	if tx > 0 or ty > 0:
		# Top hole
		var a = boundaries[0][-1]
		var b = a + num_verts
		var c = b + num_verts
		var d = c + num_verts
		connectors.append_array([a, b, c, c, d, a])
		
		# Bottom hole
		a = boundaries[2][0] + num_verts * 4
		b = a + num_verts
		c = b + num_verts
		d = c + num_verts
		connectors.append_array([a, b, c, c, d, a])

	# Side holes
	var sides = []
	if ty > 0:
		sides = [[7, 0], [1, 2], [3, 4], [5, 6]]
	
	for side in sides:
		var i = side[0]
		var j = side[1]
		
		# First side
		var patch_index = i
		var patch: int = patch_index / 2
		var next_patch = 4 + (4 - patch) % 4
		var boundary_a = boundaries[2].duplicate()
		for ii in range(boundary_a.size()):
			boundary_a[ii] += num_verts * patch
		var boundary_b = boundaries[1].duplicate()
		for ii in range(boundary_b.size()):
			boundary_b[ii] += num_verts * next_patch
		var a
		var b
		if patch_index % 2 == 0:
			a = boundary_a[0]
			b = boundary_b[n-1]
		else:
			a = boundary_a[n-1]
			b = boundary_b[0]
		
		# Second side
		patch_index = j
		patch = patch_index / 2
		next_patch = 4 + (4 - patch) % 4
		boundary_a = boundaries[2].duplicate()
		for ii in range(boundary_a.size()):
			boundary_a[ii] += num_verts * patch
		boundary_b = boundaries[1].duplicate()
		for ii in range(boundary_b.size()):
			boundary_b[ii] += num_verts * next_patch
		var c
		var d
		if patch_index % 2 == 0:
			c = boundary_a[0]
			d = boundary_b[n-1]
		else:
			c = boundary_a[n-1]
			d = boundary_b[0]
		
		connectors.append_array([a, b, d, d, c, a])

	connectors.reverse()
	return connectors

func get_connector_pairs(ndivisions: int) -> Array:
	var n = 2**ndivisions + 1
	var num_verts: int = n * (n + 1) / 2
	var tx = 1
	var ty = 1
	var tz = 1

	var boundaries: Array[PackedInt32Array] = get_boundary_indices(ndivisions)
	assert(boundaries.size() == 3)
	var connector_pairs: Array[PackedInt32Array] = []

	if radius > 0:
		# Top half
		for patch in range(4):
			if patch % 2 == 0 and tz == 0:
				continue
			if patch % 2 == 1 and tx == 0:
				continue
			var next_patch = (patch + 1) % 4
			var boundary_a = boundaries[1].duplicate()
			for i in range(boundary_a.size()):
				boundary_a[i] += num_verts * patch
			var boundary_b = boundaries[0].duplicate()
			for i in range(boundary_b.size()):
				boundary_b[i] += num_verts * next_patch
			for i in range(n-1):
				var a = boundary_a[i]
				var b = boundary_b[i]
				var c = boundary_a[i+1]
				var d = boundary_b[i+1]
				connector_pairs.append(PackedInt32Array([a, b]))
				connector_pairs.append(PackedInt32Array([c, d]))
		
		# Bottom half
		for patch in range(4, 8):
			if patch % 2 == 0 and tx == 0:
				continue
			if patch % 2 == 1 and tz == 0:
				continue
			var next_patch = 4 + (patch + 1) % 4
			var boundary_a = boundaries[0].duplicate()
			for i in range(boundary_a.size()):
				boundary_a[i] += num_verts * patch
			var boundary_b = boundaries[2].duplicate()
			for i in range(boundary_b.size()):
				boundary_b[i] += num_verts * next_patch
			for i in range(n-1):
				var a = boundary_a[i]
				var b = boundary_b[i]
				var c = boundary_a[i+1]
				var d = boundary_b[i+1]
				connector_pairs.append(PackedInt32Array([a, b]))
				connector_pairs.append(PackedInt32Array([c, d]))
		
		# Connect top patch to bottom patch
		if ty > 0:
			for patch in range(4):
				var next_patch = 4 + (4 - patch) % 4
				var boundary_a = boundaries[2].duplicate()
				for i in range(boundary_a.size()):
					boundary_a[i] += num_verts * patch
				var boundary_b = boundaries[1].duplicate()
				for i in range(boundary_b.size()):
					boundary_b[i] += num_verts * next_patch
				for i in range(n-1):
					var a = boundary_a[i]
					var b = boundary_b[n-1-i]
					var c = boundary_a[i+1]
					var d = boundary_b[n-1-i-1]
					connector_pairs.append(PackedInt32Array([a, b]))
					connector_pairs.append(PackedInt32Array([c, d]))

	if tx > 0 or ty > 0:
		# Top hole
		var a = boundaries[0][-1]
		var b = a + num_verts
		var c = b + num_verts
		var d = c + num_verts
		connector_pairs.append(PackedInt32Array([a, b, c, d]))
		
		# Bottom hole
		a = boundaries[2][0] + num_verts * 4
		b = a + num_verts
		c = b + num_verts
		d = c + num_verts
		connector_pairs.append(PackedInt32Array([a, b, c, d]))

	# Side holes
	var sides = []
	if ty > 0:
		sides = [[7, 0], [1, 2], [3, 4], [5, 6]]
	
	for side in sides:
		var i = side[0]
		var j = side[1]
		
		# First side
		var patch_index = i
		var patch: int = patch_index / 2
		var next_patch = 4 + (4 - patch) % 4
		var boundary_a = boundaries[2].duplicate()
		for ii in range(boundary_a.size()):
			boundary_a[ii] += num_verts * patch
		var boundary_b = boundaries[1].duplicate()
		for ii in range(boundary_b.size()):
			boundary_b[ii] += num_verts * next_patch
		var a
		var b
		if patch_index % 2 == 0:
			a = boundary_a[0]
			b = boundary_b[n-1]
		else:
			a = boundary_a[n-1]
			b = boundary_b[0]
		
		# Second side
		patch_index = j
		patch = patch_index / 2
		next_patch = 4 + (4 - patch) % 4
		boundary_a = boundaries[2].duplicate()
		for ii in range(boundary_a.size()):
			boundary_a[ii] += num_verts * patch
		boundary_b = boundaries[1].duplicate()
		for ii in range(boundary_b.size()):
			boundary_b[ii] += num_verts * next_patch
		var c
		var d
		if patch_index % 2 == 0:
			c = boundary_a[0]
			d = boundary_b[n-1]
		else:
			c = boundary_a[n-1]
			d = boundary_b[0]
		
		connector_pairs.append(PackedInt32Array([a, b]))
		connector_pairs.append(PackedInt32Array([c, d]))

	return connector_pairs

func get_boundary_indices(ndivisions: int):
	var n: int = 2 ** ndivisions + 1
	var boundaries: Array[PackedInt32Array] = [PackedInt32Array(),PackedInt32Array(),PackedInt32Array()]
	var a: int = 0
	var b: int = 0
	var c: int = 0
	var j0: int = 0

	# Initialize the inner arrays with the required size
	for i in range(3):
		boundaries[i].resize(n)

	var row
	for col_index in range(n - 1):
		var col_height: int = n - 1 - col_index
		var j1: int = j0 + 1

		boundaries[0][a] = j0
		a += 1

		for row_ in range(col_height - 1):
			if col_height == n - 1:
				boundaries[2][c] = j0 + row_
				c += 1

		row = col_height - 1
		if col_height == n - 1:
			boundaries[2][c] = j0 + row
			c += 1
			boundaries[2][c] = j1 + row
			c += 1

		boundaries[1][b] = j1 + row
		b += 1

		j0 += col_height + 1

	boundaries[0][a] = j0 + row
	boundaries[1][b] = j0 + row

	return boundaries
