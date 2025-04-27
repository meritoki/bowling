encapsulated package Bowling
  type Feet = Real(unit="ft");
  type Board = Real(unit="board");
  type Position = Real(unit="m");
  type Distance = Real(unit="m");
  type Velocity = Real(unit="m/s");
  type Force = Real(unit="N");
  type Torque = Real(unit = "N.m");
  type Angle = Real(unit = "rad");
  type AngularVelocity = Real(unit = "rad/s");
  type Inertia = Real(unit = "N.m");
  model Ball
      Position x(start=0.33);
      Board xBoard;
      Position y(start=0);
      Feet yFeet;
      AngularVelocity w_x(start=-9);
      AngularVelocity w_y(start=-9);
      AngularVelocity w_z(start=5);
      Velocity v_x(start=0);
      Velocity v_y(start=8);
      Velocity v_b_x;
      Velocity v_b_y;
      Velocity v_b;
      Torque t_x;
      Torque t_y;
      Torque t_z;
      Angle th(start=(1*(3.1415926535/180)));
      parameter Real RG = 6.35/100;
      parameter Real Diff = 0.1/100;
      parameter Real IntDiff = 0;
      parameter Inertia I_y = m * (RG)^2;
      parameter Inertia I_x = m * (RG + Diff)^2;
      parameter Inertia I_z = m * (RG + Diff + IntDiff)^2;
      parameter Real r = 10.85/100;
      parameter Real m = 6.8;
      Real mu(start = 0.02);
      constant Real g = 9.8;
  equation
//Equation 1
      v_b_x = v_x-r*w_y;
//Equation 2
      v_b_y = v_y+r*w_x;
//Equation 3
      v_b = sqrt(v_b_x^2+v_b_y^2);
//Implicit Equation
      der(x) = v_x;
//Implicit Equation
      der(y) = v_y;
//Equation 10
      der(v_x) = -mu*g*(v_b_x/v_b);
//Equation 11
      der(v_y) = -mu*g*(v_b_y/v_b);
//Equation 12
      der(w_x)*cos(th)+der(w_y)*sin(th) = (1/I_x)*(m*r*(der(v_y)*cos(th)+der(v_x)*sin(th))+(I_y - I_z)*(-w_x*sin(th)+w_y*cos(th))*w_z);
//Equation 13
      -der(w_x)*sin(th)+der(w_y)*cos(th) = (1/I_y)*(m*r*(-der(v_y)*sin(th)-der(v_x)*cos(th))+(I_z - I_x)*(w_x*cos(th)+w_y*sin(th))*w_z);
//Equation 14
      der(w_z) = ((I_x - I_y)/I_z)*(w_x*cos(th)+w_y*sin(th))*(-w_x*sin(th)+w_y*cos(th));
//Equation 15
      der(th) = der(atan(w_y/w_x));
  
      t_x = I_x*der(w_x)-(I_y-I_z)*w_y*w_z;
      t_y = I_y*der(w_y)-(I_z-I_x)*w_z*w_x;
      t_z = I_z*der(w_z)-(I_x-I_y)*w_x*w_y;
      
      yFeet = y * 3.28084;
      xBoard = (x * 39.3701)/1.07;
  
      if y < 12.192 then
        mu = 0.07;
      else
        mu = 0.07;
      end if;
  annotation(
      experiment(StartTime = 0, StopTime = 3, Tolerance = 1e-06, Interval = 0.01));
  end Ball;
  
end Bowling;
