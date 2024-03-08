
"""
	Position(pos_x, pos_y, pos_z)

Data structure for the position vector of one particle.

# Arguments
- `pos_x::Real`: the position of the particle on the x axis
- `pos_y::Real`: the position of the particle on the y axis
- `pos_z::Real`: the position of the particle on the z axis
"""
mutable struct Position
	x::Real
	y::Real
	z::Real
end


"""
	Velocity(vel_x, vel_y, vel_z)

Data structure for the velocity vector of one particle.

# Arguments
- `vel_x::Real`: the velocity of the particle in the x direction
- `vel_y::Real`: the velocity of the particle in the y direction
- `vel_z::Real`: the velocity of the particle in the z direction
"""
mutable struct Velocity
	x::Real
	y::Real
	z::Real
end


"""
	Acceleration(acc_x, acc_y, acc_z)

Data structure for the acceleration vector of one particle.

# Arguments
- `acc_x::Real`: the acceleration of the particle in the x direction
- `acc_y::Real`: the acceleration of the particle in the y direction
- `acc_z::Real`: the acceleration of the particle in the z direction
"""
struct Acceleration
	x::Real
	y::Real
	z::Real
end


"""
	Particle(position, velocity, acceleration, sigma, epsilon, mass)	

Data structure for particles.

# Arguments
- `position::Position`: the position vector of the particle in x, y, and z
- `velocity::Velocity`: the velocity vector of the particle in x, y, and z
- `acceleration::Acceleration`: the acceleration vector of the particle in x, y, and z 
- `sigma::Real`: the size of the particle
- `epsilon::Real`: the depth of the potential well
- `mass::Real`: the mass of the particle

"""
struct Particle 
	position::Position
	velocity::Velocity
	acceleration::Acceleration
	sigma::Real
	epsilon::Real
	mass::Real
end


"""
	Box(len_x, len_y, len_z)

Data structure for a simulation box. The box is built in the positive x, y, and z 
directions, with the origin in one corner. All positions in the main box should be
therefore positive.

# Arguments
- `len_x::Real`: The length of the box in the x direction
- `len_x::Real`: The length of the box in the y direction
- `len_x::Real`: The length of the box in the z direction
"""
struct Box
	len_x::Real
	len_y::Real
	len_z::Real
end


function lj_potential(sigma, epsilon, r)
	sigma::Real
	epsilon::Real
	r::Real

	return 4 * epsilon * ((sigma/r)^12 - (sigma/r)^6)
end


function lj_force(sigma, epsilon, r)
	sigma::Real
	epsilon::Real
	r::Real

	return -24 * epsilon/r * (2 * (sigma/r)^12 - (sigma/r)^6)
end


function generate_particles(num; sigma, epsilon, mass, box)
	num::Integer
	sigma::Real
	epsilon::Real
	mass::Real
	box::Box


	limx = box.len_x - 1.2
	limy = box.len_y - 1.2
	limz = box.len_z - 1.2

	if num * 1.728 >= limx * limy * limz
		throw(DomainError(num, "number of particles is too large for the given box"))
	end

	particles = Array{Particle}(undef, num)
	for i in range(0, num-1)
		
		posx = (i * 1.2) % limx
		j = fld(i * 1.2, limx)

		posy = (j * 1.2) % limy
		k = fld(j * 1.2, limy)

		posz = (k * 1.2) % limz

		pos = Position(posx, posy, posz)
		vel = Velocity(0, 0, 0)
		acc = Acceleration(0, 0, 0)

		particle = Particle(pos, vel, acc, sigma, epsilon, mass)
		particles[i+1] = particle
	end

	return particles
end

function main()
	box = Box(10, 10, 10)

	particles = generate_particles(10, sigma=1, epsilon=1, mass=1, box=box)
	println(particles)
end

main()