# Happy-Tree
Computer graphics group project using OpenGL

Content includes:
- Procedural tree generation: Robbie Iversen
- Wind based physics: Kieran Mckay
- Mesh to particle system conversion: Jack Purvis


# Controls

To adjust tree generation values:


To adjust wind values:

// Wind enabled by default, but inital X-axis and Z-axis values are 0
Stop/start wind: 'f'

// Wind force on the X-axis
Increase X-axis: 'j'
Decrease X-axis: 'n'

// Wind force on the Z-axis
Increase Z-axis: 'k'
Decrease Z-axis: 'm'

// This value controlls the "swing" factor of the tree
Increase wind coefficent: 'h'
Decrease wind coefficent: 'b'

// Controls how fast the tree moves by incrementing a time value by a larger or smaller amount per frame
Increase time incrementation: 'g'
Decrease time incrementation: 'v'


To adjust particle generation:

// Tree particle generation enable by defualt
Switch to bunny mode: 'e'

// Generate particles for the given mesh(s)
Generate in realtime: 'q'
Perform one iteration of the generation algorithm: 'right mouse button'
Perform 100 iterations of the generation algorithm: 'middle mouse button'
