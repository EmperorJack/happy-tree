# Happy-Tree

**Computer graphics group project using OpenGL**

- Procedural tree generation: Robbie Iversen
- Wind based physics: Kieran Mckay
- Mesh to particle system conversion: Jack Purvis

## Build Instructions

- After building the work directory with cmake
- Run `./build/bin/project` from the project directory

## Controls

**To adjust tree generation values**

All tree controls are accessed by holding: ctrl

- Increase Tree Height: 'W'
- Decrease Tree Height: 'S'
	
- Increase Branch Length: 'D'
- Decrease Branch Length: 'A'

- Increase Branch Thickness: 'E'
- Decrease Branch Thickness: 'Q'
	
- Re-generate tree: 'R'
- Toggle Line view: 'T'

**To adjust wind values**

Wind enabled by default, but inital X-axis and Z-axis values are 0

- Stop/start wind: 'f'

Wind force on the X-axis

- Increase X-axis: 'j'
- Decrease X-axis: 'n'

Wind force on the Z-axis

- Increase Z-axis: 'k'
- Decrease Z-axis: 'm'

This value controls the "swing" factor of the tree

- Increase wind coefficent: 'h'
- Decrease wind coefficent: 'b'

Controls how fast the tree moves by incrementing a time value by a larger or smaller amount per frame

- Increase time incrementation: 'g'
- Decrease time incrementation: 'v'

**To perform particle generation**

Tree particle generation enable by default

- Switch to bunny mode: 'e'

Generate particles for the given mesh(s)

- Generate in realtime: 'q'
- Perform one iteration of the generation algorithm: 'right mouse button'
- Perform 100 iterations of the generation algorithm: 'middle mouse button'

Animate the resultant particle system

- Explode the particle system once it's finished building: 'space'
