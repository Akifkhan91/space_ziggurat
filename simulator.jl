function simulator(initital_and_boundary_conditions,data,thrust_vector,Δt)
    x = zeros(Float64, data.T)
    y = zeros(Float64, data.T)
    z = zeros(Float64, data.T)

    x[1] = initital_and_boundary_conditions.initial_position[1]
    y[1] = initital_and_boundary_conditions.initial_position[2]
    z[1] = initital_and_boundary_conditions.initial_position[3]

    x[data.T] = initital_and_boundary_conditions.final_position[1]
    y[data.T] = initital_and_boundary_conditions.final_position[2]
    z[data.T] = initital_and_boundary_conditions.final_position[3]

    v_x = zeros(Float64, data.T)
    v_y = zeros(Float64, data.T)
    v_z = zeros(Float64, data.T)

    v_x[1] = initital_and_boundary_conditions.initial_velocity[1]
    v_y[1] = initital_and_boundary_conditions.initial_velocity[2]
    v_z[1] = initital_and_boundary_conditions.initial_velocity[3]

    a_x = zeros(Float64, data.T-1)
    a_y = zeros(Float64, data.T-1)
    a_z = zeros(Float64, data.T-1)
    
    for t=2:data.T
        a_x[t-1] =  -data.scaled_mu_const * 
                    x[t-1]/(x[t-1]^2 + y[t-1]^2 + z[t-1]^2)^1.5 +  
                     thrust_vector[t-1,1]/data.initial_mass
        
        a_y[t-1] = -data.scaled_mu_const * 
                    y[t-1]/(x[t-1]^2 + y[t-1]^2 + z[t-1]^2)^1.5 +  
                     thrust_vector[t-1,2]/data.initial_mass
        
        a_z[t-1] = -data.scaled_mu_const * 
                    z[t-1]/(x[t-1]^2 + y[t-1]^2 + z[t-1]^2)^1.5 +  
                     thrust_vector[t-1,3]/data.initial_mass
        
        v_x[t] = a_x[t-1]*Δt + v_x[t-1] 
        v_y[t] = a_y[t-1]*Δt + v_y[t-1]
        v_z[t] = a_z[t-1]*Δt + v_z[t-1]

        x[t] = v_x[t] * Δt + x[t-1]
        y[t] = v_y[t] * Δt + y[t-1]
        z[t] = v_z[t] * Δt + z[t-1]
        end

end
