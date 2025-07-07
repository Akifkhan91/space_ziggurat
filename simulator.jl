function simulator(initital_and_boundary_conditions,data,ord,thrust_vector,Δt)
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
    v_x[t] = v_x[t-1] - scaled_mu_constant * Δt/ (x[t-1]^2 + y[t-1]^2 + z[t-1]^2)
    v_y[t] = v_y[t-1]
    v_z[t] = v_z[t-1]

    x[t] = v_x[t] * Δt + x[t-1]
    y[t] = v_y[t] * Δt + y[t-1]
    z[t] = v_z[t] * Δt + z[t-1]

    if t == 2
        a_x[t-1] = (v_x[t] - v_x[t-1]) / Δt
        a_y[t-1] = (v_y[t] - v_y[t-1]) / Δt
        a_z[t-1] = (v_z[t] - v_z[t-1]) / Δt
    else
        a_x[t-1] = (x[t] - 2 * x[t-1] +  x[t-2]) / Δt^2
        a_y[t-1] = (y[t] - 2 * y[t-1] +  y[t-2]) / Δt^2
        a_z[t-1] = (z[t] - 2 * z[t-1] +  z[t-2]) / Δt^2
    end
    
    end

end