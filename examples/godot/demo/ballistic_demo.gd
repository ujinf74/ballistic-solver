extends Node

func _ready() -> void:
	var solver := BallisticSolver.new()
	var result := solver.solve(
		Vector3(120.0, 30.0, 5.0),
		Vector3(2.0, -1.0, 0.0),
		90.0,
		0.002
	)
	print(result)
