encapsulated package Bowling
  type Velocity = Real(unit = "v");
  type AngularVelocity = Real(unit = "w");
  type Inertia = Real(unit = "I");
  model Ball
    parameter Real r = 10.85;
    parameter Real m = 6.8;
    parameter Real RG = 6.35;
    parameter Real Diff = 0.1;
    parameter Real IntDiff = 0;
    parameter Inertia I_y = m * (RG)^2;
    parameter Inertia I_x = m * (RG + Diff)^2;
    parameter Inertia I_z = m * (RG + Diff + IntDiff)^2;
    parameter Velocity v_x = 8;
    parameter Velocity v_y = 8;
    parameter Velocity v_b_x;
    parameter Velocity v_b_y;
    parameter AngularVelocity w_x = -30;
    parameter AngularVelocity w_y = -30;
    parameter AngularVelocity w_z = 10;
  
  equation
  v_b_x = v_x - (r * w_y);
  v_b_y = v_y + (r * w_x);
  
  end Ball;
end Bowling;
