extends Node3D

# Moving-target intercept showcase. A target is tracked from NOISY, position-only
# measurements via BallisticTracker (CA-Kalman), and the launcher fires a
# continuous burst using the tracker's lead-fire solution. Rounds are simulated
# and scored against the TRUE target so green markers reflect real intercepts.
# The same machinery applies to robotics, perimeter security, agricultural
# deterrence, games/sim, and defense -- it only sees relative kinematics.
# Units are metres; the launcher sits at the origin, +Y is up.

# Godot frame: X = down-range, Y = altitude (up), Z = cross-range.
# A diving, gently curving (constant-acceleration) threat -- not a straight line --
# engaged with slower rounds so the ballistic arc / drop is clearly visible.
const TARGET_START := Vector3(480.0, 130.0, 80.0)
const TARGET_VEL := Vector3(-115.0, -16.0, -20.0)   # inbound + gentle dive + slight cross
const TARGET_ACCEL := Vector3(0.0, 0.0, 6.0)        # cross-range curve (CA, tracker-exact)
const MUZZLE_SPEED := 400.0                          # keeps time-of-flight short -> reliable lead
const K_DRAG := 0.0015
const MEAS_NOISE := 1.5        # position measurement noise std, m
const PROCESS_NOISE := 5.0     # tracker jerk spectral density (tuned for 60 Hz updates)
const FIRE_INTERVAL := 0.11
const FIRE_RANGE := 440.0
const MIN_RANGE := 90.0
const HIT_THRESH := 5.0
const RUN_TIME := 3.6
const RESET_DELAY := 2.5
const MAX_PROJECTILES := 64
const MEAS_TRAIL := 24

var solver: BallisticSolver = BallisticSolver.new()
var tracker: BallisticTracker = BallisticTracker.new()

var target: MeshInstance3D
var pred_marker: MeshInstance3D
var status_label: Label
var projectiles: Array[Dictionary] = []
var meas_markers: Array[MeshInstance3D] = []
var hit_markers: Array[MeshInstance3D] = []
var est_marker: MeshInstance3D
var lead_line: MeshInstance3D

var elapsed := 0.0
var next_fire := 0.0
var fired := 0
var hits := 0
var best_miss := INF


func true_pos(t: float) -> Vector3:
	return TARGET_START + TARGET_VEL * t + 0.5 * TARGET_ACCEL * t * t


func _ready() -> void:
	_build_scene()
	_reset()


func _process(delta: float) -> void:
	elapsed += delta
	var tp := true_pos(elapsed)
	target.position = tp

	var meas := tp + Vector3(randfn(0.0, MEAS_NOISE), randfn(0.0, MEAS_NOISE), randfn(0.0, MEAS_NOISE))
	tracker.update(elapsed, meas)
	_spawn_meas_marker(meas)

	var rng := tp.length()
	if rng < FIRE_RANGE and rng > MIN_RANGE and elapsed >= next_fire and elapsed <= RUN_TIME:
		_fire()
		next_fire = elapsed + FIRE_INTERVAL

	# tracking story: noisy measurement (red) -> denoised estimate (cyan) -> lead (gold)
	est_marker.position = tracker.estimated_position()
	var lim := lead_line.mesh as ImmediateMesh
	lim.clear_surfaces()
	if pred_marker.position != Vector3.ZERO:
		lim.surface_begin(Mesh.PRIMITIVE_LINES)
		lim.surface_add_vertex(Vector3.ZERO)
		lim.surface_add_vertex(pred_marker.position)
		lim.surface_end()

	_update_projectiles()
	_update_status(rng)

	if elapsed > RUN_TIME + RESET_DELAY:
		_reset()


func _fire() -> void:
	var sol: Dictionary = Dictionary(tracker.solve(MUZZLE_SPEED, K_DRAG))
	if not bool(sol.get("success", false)):
		return

	var theta := float(sol["theta"])
	var phi := float(sol["phi"])
	var t_star := float(sol["t_star"])
	pred_marker.position = tracker.predict(t_star)

	var steps: int = max(12, int(ceil(t_star / 0.02)))
	var pts := PackedVector3Array(solver.simulate_from_angles(
		MUZZLE_SPEED, theta, phi, K_DRAG, t_star, t_star / float(steps)))
	if pts.is_empty():
		return

	var node := _sphere(2.2, Color(0.6, 0.7, 0.8))
	add_child(node)
	var trail := MeshInstance3D.new()
	trail.mesh = ImmediateMesh.new()
	var tmat := StandardMaterial3D.new()
	tmat.shading_mode = BaseMaterial3D.SHADING_MODE_UNSHADED
	tmat.albedo_color = Color(1.0, 0.72, 0.26)
	tmat.emission_enabled = true
	tmat.emission = Color(1.0, 0.72, 0.26)
	trail.material_override = tmat
	add_child(trail)
	projectiles.append({
		"node": node,
		"trail": trail,
		"points": pts,
		"start_time": elapsed,
		"t_star": t_star,
		"min_dist": INF,
		"hit": false,
	})
	fired += 1
	if projectiles.size() > MAX_PROJECTILES:
		_remove_projectile(0)


func _update_projectiles() -> void:
	for i in range(projectiles.size() - 1, -1, -1):
		var info: Dictionary = projectiles[i]
		var node: MeshInstance3D = info["node"] as MeshInstance3D
		var pts := PackedVector3Array(info["points"])
		var t_star := float(info["t_star"])
		var age := elapsed - float(info["start_time"])

		if age > t_star + 0.35:
			_remove_projectile(i)
			continue

		var pos := _sample_points(pts, age, t_star)
		node.position = pos
		var im := (info["trail"] as MeshInstance3D).mesh as ImmediateMesh
		var f := clampf(age / maxf(0.001, t_star), 0.0, 1.0) * float(pts.size() - 1)
		var upto := int(floor(f))
		im.clear_surfaces()
		if upto >= 1:
			im.surface_begin(Mesh.PRIMITIVE_LINE_STRIP)
			for k in range(upto + 1):
				im.surface_add_vertex(pts[k])
			im.surface_add_vertex(pos)
			im.surface_end()
		var d := pos.distance_to(true_pos(elapsed))
		if d < float(info["min_dist"]):
			info["min_dist"] = d
		if not bool(info["hit"]) and d <= HIT_THRESH:
			info["hit"] = true
			hits += 1
			best_miss = minf(best_miss, d)
			_spawn_hit_marker(pos)
			_set_color(node, Color(0.2, 1.0, 0.3))


func _update_status(rng: float) -> void:
	var pk := 0.0
	if fired > 0:
		pk = 100.0 * float(hits) / float(fired)
	status_label.text = "moving-target intercept  t=%.1fs  range=%.0fm  fired=%d  hits<%.0fm=%d (%.0f%%)  best=%.2fm" % [
		elapsed, rng, fired, HIT_THRESH, hits, pk,
		(best_miss if best_miss < INF else 0.0)
	]


func _reset() -> void:
	for i in range(projectiles.size() - 1, -1, -1):
		_remove_projectile(i)
	for m in meas_markers:
		m.queue_free()
	meas_markers.clear()
	for m in hit_markers:
		m.queue_free()
	hit_markers.clear()
	tracker.configure(PROCESS_NOISE, MEAS_NOISE * MEAS_NOISE)
	elapsed = 0.0
	next_fire = 0.0
	fired = 0
	hits = 0
	best_miss = INF
	target.position = true_pos(0.0)
	pred_marker.position = Vector3.ZERO


func _spawn_meas_marker(p: Vector3) -> void:
	var m := _sphere(1.6, Color(1.0, 0.54, 0.43))
	m.position = p
	add_child(m)
	meas_markers.append(m)
	if meas_markers.size() > MEAS_TRAIL:
		var old: MeshInstance3D = meas_markers.pop_front()
		old.queue_free()


func _spawn_hit_marker(p: Vector3) -> void:
	var m := _sphere(2.6, Color(0.2, 1.0, 0.3))
	m.position = p
	add_child(m)
	hit_markers.append(m)
	if hit_markers.size() > 64:
		var old: MeshInstance3D = hit_markers.pop_front()
		old.queue_free()


func _remove_projectile(index: int) -> void:
	var info: Dictionary = projectiles[index]
	(info["node"] as MeshInstance3D).queue_free()
	(info["trail"] as MeshInstance3D).queue_free()
	projectiles.remove_at(index)


func _build_scene() -> void:
	var camera := Camera3D.new()
	camera.fov = 50.0
	camera.fov = 58.0
	camera.look_at_from_position(Vector3(-170.0, 135.0, 205.0), Vector3(210.0, 45.0, 45.0), Vector3.UP)
	camera.current = true
	add_child(camera)

	var light := DirectionalLight3D.new()
	light.rotation_degrees = Vector3(-55.0, 35.0, 0.0)
	add_child(light)

	var ground := MeshInstance3D.new()
	var ground_mesh := PlaneMesh.new()
	ground_mesh.size = Vector2(1100.0, 700.0)
	ground.mesh = ground_mesh
	ground.position = Vector3(310.0, 0.0, 60.0)
	var ground_mat := StandardMaterial3D.new()
	ground_mat.albedo_color = Color(0.10, 0.13, 0.16)
	ground.material_override = ground_mat
	add_child(ground)

	var gun := _sphere(6.0, Color(0.95, 0.95, 0.95))
	gun.position = Vector3.ZERO
	add_child(gun)

	target = _sphere(5.0, Color(1.0, 0.25, 0.18))
	add_child(target)

	pred_marker = _sphere(4.0, Color(1.0, 0.82, 0.16))
	add_child(pred_marker)

	# denoised track estimate (the CA-Kalman position) -- cyan
	est_marker = _sphere(3.2, Color(0.3, 0.8, 1.0))
	add_child(est_marker)

	# lead line: launcher -> predicted intercept (the firing solution), gold
	lead_line = MeshInstance3D.new()
	lead_line.mesh = ImmediateMesh.new()
	var lmat := StandardMaterial3D.new()
	lmat.shading_mode = BaseMaterial3D.SHADING_MODE_UNSHADED
	lmat.albedo_color = Color(1.0, 0.82, 0.16)
	lmat.emission_enabled = true
	lmat.emission = Color(1.0, 0.82, 0.16)
	lead_line.material_override = lmat
	add_child(lead_line)

	# radar-style range rings on the ground: scale + depth + theme
	var rings := MeshInstance3D.new()
	var rim := ImmediateMesh.new()
	rings.mesh = rim
	var rmat := StandardMaterial3D.new()
	rmat.shading_mode = BaseMaterial3D.SHADING_MODE_UNSHADED
	rmat.transparency = BaseMaterial3D.TRANSPARENCY_ALPHA
	rmat.albedo_color = Color(0.32, 0.5, 0.62, 0.45)
	rings.material_override = rmat
	add_child(rings)
	rim.surface_begin(Mesh.PRIMITIVE_LINES)
	for ring_r in [100, 200, 300, 400]:
		var prev := Vector3(ring_r, 0.5, 0.0)
		for k in range(1, 73):
			var a := TAU * float(k) / 72.0
			var p := Vector3(ring_r * cos(a), 0.5, ring_r * sin(a))
			rim.surface_add_vertex(prev)
			rim.surface_add_vertex(p)
			prev = p
	rim.surface_end()

	var canvas := CanvasLayer.new()
	add_child(canvas)
	status_label = Label.new()
	status_label.position = Vector2(18.0, 14.0)
	status_label.add_theme_font_size_override("font_size", 17)
	canvas.add_child(status_label)

	var legend := RichTextLabel.new()
	legend.bbcode_enabled = true
	legend.fit_content = true
	legend.scroll_active = false
	legend.position = Vector2(18.0, 44.0)
	legend.custom_minimum_size = Vector2(560.0, 0.0)
	legend.add_theme_font_size_override("normal_font_size", 15)
	legend.text = "[color=#ff5a40]● threat[/color]   [color=#ff8a6e]● noisy radar[/color]   [color=#4ccfff]● track estimate[/color]   [color=#ffd23f]● lead[/color]   [color=#ffb347]— round[/color]   [color=#33ff55]★ hit[/color]"
	canvas.add_child(legend)


func _sample_points(points: PackedVector3Array, t: float, t_star: float) -> Vector3:
	if points.is_empty():
		return Vector3.ZERO
	var f := clampf(t / maxf(0.001, t_star), 0.0, 1.0) * float(points.size() - 1)
	var i := int(floor(f))
	var j := mini(i + 1, points.size() - 1)
	return points[i].lerp(points[j], f - float(i))


func _sphere(radius: float, color: Color) -> MeshInstance3D:
	var node := MeshInstance3D.new()
	var mesh := SphereMesh.new()
	mesh.radius = radius
	mesh.height = radius * 2.0
	node.mesh = mesh
	var mat := StandardMaterial3D.new()
	mat.albedo_color = color
	mat.emission_enabled = true
	mat.emission = color * 0.3
	node.material_override = mat
	return node


func _set_color(node: MeshInstance3D, color: Color) -> void:
	var mat := StandardMaterial3D.new()
	mat.albedo_color = color
	mat.emission_enabled = true
	mat.emission = color * 0.4
	node.material_override = mat
