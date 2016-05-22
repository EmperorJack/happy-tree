# Adds the resources as targets so they show up in IDEs

# TODO List your shader source files here
SET(SHADERS
	"shaders/gouraudShader.vert"
	"shaders/gouraudShader.frag"
	"shaders/phongShader.vert"
	"shaders/phongShader.frag"
)

add_custom_target(
	res
	SOURCES
	${SHADERS}
)

source_group("shader" FILES ${SHADERS})
