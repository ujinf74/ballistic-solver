extends Node3D

const TARGET_START := Vector3(120.0, 10.0, -30.0)
const TARGET_VEL := Vector3(1.0, 0.0, 4.0)
const MUZZLE_SPEED := 120.0
const K_DRAG := 0.01
const TIME_STEP := 0.04
const FIRE_INTERVAL := 0.28
const FIRE_DURATION := 6.0
const RESET_DELAY := 2.0
const MAX_PROJECTILES := 32

var solver: BallisticSolver = BallisticSolver.new()
var last_result: Dictionary = {}
var trajectory_line: MeshInstance3D
var lead_line: MeshInstance3D
var target: MeshInstance3D
var impact_marker: MeshInstance3D
var status_label: Label
var projectiles: Array[Dictionary] = []
var elapsed_time := 0.0
var next_fire_time := 0.0
var max_flight_time := 0.0
var shots_fired := 0

func _ready() -> void:
	_build_scene()
	_reset_demo()

func _process(delta: float) -> void:
	elapsed_time += delta
	target.position = TARGET_START + TARGET_VEL * elapsed_time
	_fire_due_projectiles()
	_update_projectiles()
	_update_status()

	if elapsed_time > FIRE_DURATION + maxf(max_flight_time, 1.0) + RESET_DELAY:
		_reset_demo()

func _build_scene() -> void:
	var camera := Camera3D.new()
	camera.fov = 42.0
	camera.look_at_from_position(Vector3(40.0, 10.0, 100.0), Vector3(60.0, 0.0, 0.0), Vector3.UP)
	camera.current = true
	add_child(camera)

	var light := DirectionalLight3D.new()
	light.rotation_degrees = Vector3(-55.0, 35.0, 0.0)
	add_child(light)

	var ground := MeshInstance3D.new()
	var ground_mesh := PlaneMesh.new()
	ground_mesh.size = Vector2(180.0, 140.0)
	ground.mesh = ground_mesh
	ground.position = Vector3(65.0, -0.05, 0.0)
	var ground_mat := StandardMaterial3D.new()
	ground_mat.albedo_color = Color(0.12, 0.16, 0.18)
	ground.material_override = ground_mat
	add_child(ground)

	trajectory_line = MeshInstance3D.new()
	add_child(trajectory_line)

	lead_line = MeshInstance3D.new()
	add_child(lead_line)

	target = _sphere(1.0, Color(1.0, 0.25, 0.18))
	add_child(target)

	impact_marker = _sphere(0.5, Color(0.2, 1.0, 0.38))
	add_child(impact_marker)

	var muzzle := _sphere(1.0, Color(0.95, 0.95, 0.95))
	muzzle.position = Vector3.ZERO
	add_child(muzzle)

	var canvas := CanvasLayer.new()
	add_child(canvas)
	status_label = Label.new()
	status_label.position = Vector2(18.0, 16.0)
	status_label.add_theme_font_size_override("font_size", 18)
	canvas.add_child(status_label)

func _reset_demo() -> void:
	_clear_projectiles()
	elapsed_time = 0.0
	next_fire_time = 0.0
	max_flight_time = 0.0
	shots_fired = 0
	last_result = {}
	target.position = TARGET_START
	impact_marker.position = TARGET_START
	trajectory_line.mesh = null
	lead_line.mesh = null
	_fire_due_projectiles()
	_update_status()

func _fire_due_projectiles() -> void:
	var fire_until: float = minf(elapsed_time, FIRE_DURATION)
	while next_fire_time <= fire_until:
		_spawn_projectile(next_fire_time)
		next_fire_time += FIRE_INTERVAL

func _spawn_projectile(fire_time: float) -> void:
	var target_at_fire: Vector3 = TARGET_START + TARGET_VEL * fire_time
	var shot_result: Dictionary = Dictionary(solver.solve(target_at_fire, TARGET_VEL, MUZZLE_SPEED, K_DRAG))
	last_result = shot_result
	print(shot_result)

	if not bool(shot_result.get("success", false)):
		return

	var theta: float = float(shot_result["theta"])
	var phi: float = float(shot_result["phi"])
	var t_star: float = float(shot_result["t_star"])
	var impact: Vector3 = target_at_fire + TARGET_VEL * t_star
	impact_marker.position = impact
	max_flight_time = maxf(max_flight_time, t_star)

	var steps: int = max(8, int(ceil(t_star / TIME_STEP)))
	var points: PackedVector3Array = PackedVector3Array(solver.simulate_from_angles(
		MUZZLE_SPEED,
		theta,
		phi,
		K_DRAG,
		t_star,
		t_star / float(steps)
	))

	var projectile_node: MeshInstance3D = _sphere(0.5, Color(0.15, 0.55, 1.0))
	add_child(projectile_node)

	var line_node := MeshInstance3D.new()
	add_child(line_node)
	_set_line_mesh(line_node, points, Color(0.1, 0.85, 1.0, 0.38))

	_set_line_mesh(trajectory_line, points, Color(0.1, 0.85, 1.0))
	_set_line_mesh(lead_line, PackedVector3Array([Vector3.ZERO, impact]), Color(1.0, 0.82, 0.16))

	projectiles.append({
		"node": projectile_node,
		"line": line_node,
		"start_time": fire_time,
		"points": points,
		"t_star": t_star
	})
	shots_fired += 1

	if projectiles.size() > MAX_PROJECTILES:
		_remove_projectile(0)

func _update_projectiles() -> void:
	for i in range(projectiles.size() - 1, -1, -1):
		var info: Dictionary = projectiles[i]
		var projectile_node: MeshInstance3D = info["node"] as MeshInstance3D
		var points: PackedVector3Array = PackedVector3Array(info["points"])
		var t_star: float = float(info["t_star"])
		var age: float = elapsed_time - float(info["start_time"])

		if age > t_star + 0.45:
			_remove_projectile(i)
		else:
			projectile_node.position = _sample_points(points, age, t_star)

func _update_status() -> void:
	if last_result.is_empty():
		status_label.text = "firing %.1f s at %.2f s interval" % [FIRE_DURATION, FIRE_INTERVAL]
		return

	if not bool(last_result.get("success", false)):
		status_label.text = "solve failed: status=%s rc=%s %s" % [
			last_result.get("status", "?"),
			last_result.get("call_rc", "?"),
			last_result.get("message", "")
		]
		return

	status_label.text = "burst %.1f/%.1f s  shots=%d  active=%d  miss=%.4f m  t*=%.3f s  iter=%d" % [
		minf(elapsed_time, FIRE_DURATION),
		FIRE_DURATION,
		shots_fired,
		projectiles.size(),
		float(last_result["miss"]),
		float(last_result["t_star"]),
		int(last_result["iterations"])
	]

func _sample_points(points: PackedVector3Array, t: float, t_star: float) -> Vector3:
	if points.is_empty():
		return Vector3.ZERO

	var f: float = clampf(t / maxf(0.001, t_star), 0.0, 1.0) * float(points.size() - 1)
	var i: int = int(floor(f))
	var j: int = mini(i + 1, points.size() - 1)

	return points[i].lerp(points[j], f - float(i))

func _remove_projectile(index: int) -> void:
	var info: Dictionary = projectiles[index]
	var projectile_node: MeshInstance3D = info["node"] as MeshInstance3D
	var line_node: MeshInstance3D = info["line"] as MeshInstance3D
	projectile_node.queue_free()
	line_node.queue_free()
	projectiles.remove_at(index)

func _clear_projectiles() -> void:
	for i in range(projectiles.size() - 1, -1, -1):
		_remove_projectile(i)

func _sphere(radius: float, color: Color) -> MeshInstance3D:
	var node := MeshInstance3D.new()
	var mesh := SphereMesh.new()
	mesh.radius = radius
	mesh.height = radius * 2.0
	node.mesh = mesh
	var mat := StandardMaterial3D.new()
	mat.albedo_color = color
	mat.emission_enabled = true
	mat.emission = color * 0.25
	node.material_override = mat
	return node

func _set_line_mesh(node: MeshInstance3D, points: PackedVector3Array, color: Color) -> void:
	var mesh := ImmediateMesh.new()
	var mat := StandardMaterial3D.new()
	mat.shading_mode = BaseMaterial3D.SHADING_MODE_UNSHADED
	mat.albedo_color = color
	if color.a < 1.0:
		mat.transparency = BaseMaterial3D.TRANSPARENCY_ALPHA

	mesh.surface_begin(Mesh.PRIMITIVE_LINE_STRIP, mat)
	for p in points:
		mesh.surface_add_vertex(p)
	mesh.surface_end()
	node.mesh = mesh
