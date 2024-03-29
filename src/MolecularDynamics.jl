
"""
    Position(x, y, z)

Data structure for the position vector of one particle.

# Arguments
- `x::Real`: the position of the particle on the x axis
- `y::Real`: the position of the particle on the y axis
- `z::Real`: the position of the particle on the z axis
"""
mutable struct Position
    x::Real
    y::Real
    z::Real
end


"""
    Velocity(x, y, z)

Data structure for the velocity vector of one particle.

# Arguments
- `x::Real`: the velocity of the particle in the x direction
- `y::Real`: the velocity of the particle in the y direction
- `z::Real`: the velocity of the particle in the z direction
"""
mutable struct Velocity
    x::Real
    y::Real
    z::Real
end


"""
    Acceleration(x, y, z)

Data structure for the acceleration vector of one particle.

# Arguments
- `x::Real`: the acceleration of the particle in the x direction
- `y::Real`: the acceleration of the particle in the y direction
- `z::Real`: the acceleration of the particle in the z direction
"""
mutable struct Acceleration
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
    lenx::Real
    leny::Real
    lenz::Real
end


""" 
    generate_particles(num; sigma, epsilon, mass, box)

Generate particles inside a box and assign them initial conditions. The particles are 
generated in an ordered fashion, first filling the x-axis, then the y-axis and then the
z-axis. If the number of particles is too large for the simulation box throw and error.
"""
function generate_particles(num; sigma, epsilon, mass, box)
	num::Integer
	sigma::Real
	epsilon::Real
	mass::Real
	box::Box

    # Compute the effective size of the box
    effsize = 1.2 * sigma
	limx = box.lenx - effsize
	limy = box.leny - effsize
	limz = box.lenz - effsize

    # Compute the number of particles that fit per dimension
    spotsx = fld(limx, effsize)
    spotsy = fld(limy, effsize)
    spotsz = fld(limz, effsize)

    # Check if the number of particles fit inside the volume
    if num > spotsx * spotsy * spotsz
        err = "Number of particles used is too large for the box size!"
		throw(DomainError(num, err))
	end
    
    # Allocate particles
	particles = Array{Particle}(undef, num)
    for i in 0:1:(num-1)
		
        posx = i % spotsx * effsize
        j = fld(i, spotsx)

		posy = j % spotsy * effsize
        k = fld(j, spotsy)

        posz = k % spotsz * effsize

		pos = Position(posx, posy, posz)
		vel = Velocity(0.0, 0.0, 0.0)
		acc = Acceleration(0.0, 0.0, 0.0)

		particle = Particle(pos, vel, acc, sigma, epsilon, mass)
		particles[i+1] = particle
	end

	return particles
end


function lj_potential(sigma, epsilon, distance, cutoff)
	sigma::Real
	epsilon::Real
	distance::Real
    cutoff::Real    

    if cutoff != 0.0 && distance >= cutoff
        return 0.0
    else
        fraction = (sigma / distance)^6
        lj = 4 * epsilon * (fraction^2 - fraction)

        if cutoff == 0.0
	        return lj
        elseif distance < cutoff
            cutfraction = (sigma / cutoff)^6
            cutlj = 4 * epsilon * (cutfraction^2 - cutfraction)
            return lj - cutlj
        end
    end
end


function lj_force(sigma, epsilon, distance, dimension, cutoff)
	sigma::Real
	epsilon::Real
	distance::Real
    dimension::Real
    cutoff::Real
    
    if cutoff != 0.0 && distance >= cutoff
        return 0.0
    else 
        fraction = sigma / distance
        fraction7 = fraction^7   
	    return 24 * epsilon * dimension * (2 * fraction7^2 - fraction7 * fraction)
    end
end


function verlet(parts, oldparts, timestep, box)
    parts::Array{Particle}
    oldparts::Array{Particle}
    timestep::Real
    box::Box

    num = length(parts)
    newparts = deepcopy(parts)
    dt2 = timestep^2
    itdt = 0.5 / timestep
    
    # Keep track of energy conversion
    kinetic = 0.0
    potential = 0.0

    # Loop over all particle pairs
    for i in 1:1:num
        
        forcex = 0.0
        forcey = 0.0
        forcez = 0.0
        for j in 1:1:num    
            if i != j

                # Compute distance between particles
                dx = parts[i].position.x - parts[j].position.x
                dy = parts[i].position.y - parts[j].position.y
                dz = parts[i].position.z - parts[j].position.z

                # Compute the nearest image
                dx = dx - box.lenx * round(dx / box.lenx)
                dy = dy - box.leny * round(dy / box.leny)
                dz = dz - box.lenz * round(dz / box.lenz)
                
                dist = sqrt(dx^2 + dy^2 + dz^2)

                # Use mixing rules
                sigma = (parts[i].sigma + parts[j].sigma) * 0.5
                epsilon = sqrt(parts[i].epsilon * parts[j].epsilon)
                
                # Compute the forces in each direction
                forcex += lj_force(sigma, epsilon, dist, dx)
                forcey += lj_force(sigma, epsilon, dist, dy)
                forcez += lj_force(sigma, epsilon, dist, dz)

                # Compute potential energy
                potential += lj_potential(sigma, epsilon, dist) * 0.5
            end
        end
        
        # Compute the new acceleartions
        newparts[i].acceleration.x = forcex / parts[i].mass
        newparts[i].acceleration.y = forcey / parts[i].mass
        newparts[i].acceleration.z = forcez / parts[i].mass
        
        # Compute updated positions at t + dt
        newparts[i].position.x = 2.0 * parts[i].position.x - oldparts[i].position.x + 
                                newparts[i].acceleration.x * dt2
        
        newparts[i].position.y = 2.0 * parts[i].position.y - oldparts[i].position.y + 
                                newparts[i].acceleration.y * dt2
        
        newparts[i].position.z = 2.0 * parts[i].position.z - oldparts[i].position.z + 
                                newparts[i].acceleration.z * dt2

        # Compute the new velocities
        newparts[i].velocity.x = (newparts[i].position.x - oldparts[i].position.x) * itdt
        newparts[i].velocity.y = (newparts[i].position.y - oldparts[i].position.y) * itdt
        newparts[i].velocity.z = (newparts[i].position.z - oldparts[i].position.z) * itdt
        
        # Compute total kinetic energy
        kinetic += (newparts[i].velocity.x^2 + newparts[i].velocity.y^2 +
            newparts[i].velocity.z^2) * newparts[i].mass * 0.5
    end

    # Periodic boundary conditions
    for i in 1:1:num
        offsetx = box.lenx * floor(newparts[i].position.x / box.lenx)
        offsety = box.leny * floor(newparts[i].position.y / box.leny)
        offsetz = box.lenz * floor(newparts[i].position.z / box.lenz)
        
        parts[i].position.x = parts[i].position.x - offsetx 
        parts[i].position.y = parts[i].position.y - offsety
        parts[i].position.z = parts[i].position.z - offsetz
        
        newparts[i].position.x = newparts[i].position.x - offsetx
        newparts[i].position.y = newparts[i].position.y - offsety 
        newparts[i].position.z = newparts[i].position.z - offsetz
    end
    
    println("Potential: $potential, Kinetic: $kinetic, Total: $(potential+kinetic)")

    return newparts, parts
end


function velocity_verlet(parts, timestep, box, periodic, cutoff)
    parts::Array{Particle}
    timestep::Real
    box::Box
    periodic::Bool
    cutoff::Real

    num = length(parts)
    newparts = deepcopy(parts)
    dt2 = timestep^2

    # Keep track of energy conversion
    kinetic = 0.0
    potential = 0.0

    # Update new positions
    for i in 1:1:num
        
        # Update the position for t + dt
        newparts[i].position.x = parts[i].position.x + parts[i].velocity.x * timestep +
            0.5 * parts[i].acceleration.x * dt2

        newparts[i].position.y = parts[i].position.y + parts[i].velocity.y * timestep +
            0.5 * parts[i].acceleration.y * dt2
        
        newparts[i].position.z = parts[i].position.z + parts[i].velocity.z * timestep +
            0.5 * parts[i].acceleration.z * dt2
    end

    # Loop over all particle pairs
    for i in 1:1:num
        forcex = 0.0
        forcey = 0.0
        forcez = 0.0
        for j in 1:1:num    
            if i != j

                # Compute distance between particles
                dx = newparts[i].position.x - newparts[j].position.x
                dy = newparts[i].position.y - newparts[j].position.y
                dz = newparts[i].position.z - newparts[j].position.z

                # Compute the nearest image
                if periodic == true
                    dx = dx - box.lenx * round(dx / box.lenx)
                    dy = dy - box.leny * round(dy / box.leny)
                    dz = dz - box.lenz * round(dz / box.lenz)
                end
                
                dist = sqrt(dx^2 + dy^2 + dz^2)
                    
                # Use mixing rules
                sigma = (newparts[i].sigma + newparts[j].sigma) * 0.5
                epsilon = sqrt(newparts[i].epsilon * newparts[j].epsilon)
                
                # Compute the forces in each direction
                forcex += lj_force(sigma, epsilon, dist, dx, cutoff)
                forcey += lj_force(sigma, epsilon, dist, dy, cutoff)
                forcez += lj_force(sigma, epsilon, dist, dz, cutoff)

                # Compute total potential energy
                potential += lj_potential(sigma, epsilon, dist, cutoff) * 0.5
            end
        end
        
        # Compute the new acceleartions
        newparts[i].acceleration.x = forcex / newparts[i].mass
        newparts[i].acceleration.y = forcey / newparts[i].mass
        newparts[i].acceleration.z = forcez / newparts[i].mass
        
        # Compute the new velocities
        newparts[i].velocity.x = parts[i].velocity.x + 0.5 * (parts[i].acceleration.x +
            newparts[i].acceleration.x) * timestep

        newparts[i].velocity.y = parts[i].velocity.y + 0.5 * (parts[i].acceleration.y +
            newparts[i].acceleration.y) * timestep

        newparts[i].velocity.z = parts[i].velocity.z + 0.5 * (parts[i].acceleration.z +
            newparts[i].acceleration.z) * timestep
        
        # Compute total kinetic energy
        kinetic += (newparts[i].velocity.x^2 + newparts[i].velocity.y^2 +
            newparts[i].velocity.z^2) * newparts[i].mass * 0.5
    end
    
    # Periodic boundary conditions
    if periodic == true
        for i in 1:1:num
            offsetx = box.lenx * floor(newparts[i].position.x / box.lenx)
            offsety = box.leny * floor(newparts[i].position.y / box.leny)
            offsetz = box.lenz * floor(newparts[i].position.z / box.lenz)
        
            parts[i].position.x = parts[i].position.x - offsetx 
            parts[i].position.y = parts[i].position.y - offsety
            parts[i].position.z = parts[i].position.z - offsetz
        
            newparts[i].position.x = newparts[i].position.x - offsetx
            newparts[i].position.y = newparts[i].position.y - offsety 
            newparts[i].position.z = newparts[i].position.z - offsetz
        end
    end
    
    println("Potential: $potential, Kinetic: $kinetic, Total: $(potential+kinetic)")

    return newparts
end



function main()
	box = Box(6.0, 6.0, 6.0)
    timestep = 0.01
    start = 0
    stop = 4 
    num = 61
    cutoff = 1.4

	particles = generate_particles(num, sigma=1.0, epsilon=4.0, mass=1.0, box=box)
    
    file = open("results.xyz", "w")
    
    newparticles = deepcopy(particles)
    @time for _ in start:timestep:stop
        
        # Apply Verlet integration
        # newparticles, particles = verlet(newparticles, particles, timestep, box)
        
        # Apply Velocity Verlet integrator
        newparticles = velocity_verlet(newparticles, timestep, box, true, cutoff)
        
        # Print the results in .xyz format
        write(file, "$num \n")
        write(file, "\n")
        for (_, val) in enumerate(newparticles)
            posx = Float32(val.position.x)
            posy = Float32(val.position.y)
            posz = Float32(val.position.z)
            write(file, "P $(posx) $(posy) $(posz) \n")
        end
    end

    close(file)
end

main()
