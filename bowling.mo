encapsulated package Bowling
  type Position = Real(unit="m");
  type Distance = Real(unit="m");
  type Velocity = Real(unit="m/s");
  type Force = Real(unit="N");
  type Torque = Real(unit = "N.m");
  type Angle = Real(unit = "rad");
  type AngularVelocity = Real(unit = "rad/s");
  type Inertia = Real(unit = "N.m");
  model Ball
      Position x(start=0.83);
      Position y(start=0);
      AngularVelocity w_x(start=-30);
      AngularVelocity w_y(start=-30);
      AngularVelocity w_z(start=10);
      Velocity v_x(start=0);
      Velocity v_y(start=8);
      Velocity v_b_x;
      Velocity v_b_y;
      Velocity v_b;
      Angle th(start=(1.8*(3.1415926535/180)));
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
  
      v_b_x = v_x-r*w_y;
      v_b_y = v_y+r*w_x;
      v_b = sqrt(v_b_x^2+v_b_y^2);
      der(x) = v_x;
      der(y) = v_y;
      der(v_x) = -mu*g*(v_b_x/v_b);
      der(v_y) = -mu*g*(v_b_y/v_b);
      //der(th)=(der(w_y)*w_x-der(w_x)*w_y)/(w_y^2+w_x^2);
      der(w_x)*cos(th)+der(w_y)*sin(th) = (1/I_x)*(m*r*(der(v_y)*cos(th)+der(v_x)*sin(th))+(I_y - I_z)*(-w_x*sin(th)+w_y*cos(th))*w_z);
      //Equation 13
      -der(w_x)*sin(th)+der(w_y)*cos(th) = (1/I_y)*(m*r*(-der(v_y)*sin(th)-der(v_x)*cos(th))+(I_z - I_x)*(w_x*cos(th)+w_y*sin(th))*w_z);
      //Equation 14
      der(w_z) = ((I_x - I_y)/I_z)*(w_x*cos(th)+w_y*sin(th))*(-w_x*sin(th)+w_y*cos(th));
      der(th) = der(atan(w_y/w_x));
      //der(th)=(der(w_y)*w_x-der(w_x)*w_y)/(w_y^2+w_x^2);
  
          if y < 12.192 then
        mu = 0.04;
      else
        mu = 0.2;
      end if;
  annotation(
      experiment(StartTime = 0, StopTime = 5, Tolerance = 1e-06, Interval = 0.01));
  end Ball;
  
end Bowling;
