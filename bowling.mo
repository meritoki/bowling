encapsulated package Bowling
  type Feet = Real(unit="ft");
  type Board = Real(unit="board");
  type Position = Real(unit="m");
  type Distance = Real(unit="m");
  type Velocity = Real(unit="m/s");
  type MPH = Real(unit="mph");
  type Force = Real(unit="N");
  type Torque = Real(unit = "N.m");
  type Angle = Real(unit = "rad");
  type Deg = Real(unit = "deg");
  type AngularVelocity = Real(unit = "rad/s");
  type Inertia = Real(unit = "N.m");
  model Ball
      //constant Real launch_angle = 3.0;
      //constant Angle th = launch_angle*(3.1415926535/180);
      //Velocity ball_speed = 8.0;
      Position x(start=((39-xBoard0)/36.794)); // 0.53 is 19.5 (center of 20 board)
      //Position x(start=0.53); // 0.53 is 20 board
      Board xBoard;
      Position y(start=0);
      Feet yFeet;
      AngularVelocity w_x(start=-30);
      AngularVelocity w_y(start=-30);
      AngularVelocity w_z(start=10);
      Angle phi(start=45*(3.1415926535/180)); // 65.9 (optimal is 33.85 deg)
      Velocity v_x(start=Ball_Speed*(0.44704)*sin(Launch_Angle*(3.1415926535/180))); // 1 mph = 0.44704 m/s
      Velocity v_y(start=Ball_Speed*(0.44704)*cos(Launch_Angle*(3.1415926535/180)));
      //Velocity v_x(start=0.4187); // launch 3 deg
      //Velocity v_y(start=7.989);  // launch 3 Deg
      //Velocity v_x(start=0.5598);   // launch 4 deg
      //Velocity v_y(start=7.9804);   // launch 4 Deg
      //Velocity v_x(start=0.8388);   // launch 6 deg
      //Velocity v_y(start=7.9559);   // launch 6 Deg
      Velocity v_b_x;
      Velocity v_b_y;
      Velocity v_b;    
      parameter Board xBoard0 = 19.5;  // Board ball crosses at foul line
      parameter Deg Launch_Angle = 3;
      parameter MPH Ball_Speed = 17.9; // (8 m/s)
      parameter Real RG = 6.35/100;
      parameter Real Diff = 0.1/100;
      parameter Real IntDiff = 0;
      parameter Inertia I_y_p = m * (RG)^2;
      parameter Inertia I_x_p = m * (RG + Diff)^2;
      parameter Inertia I_z_p = m * (RG + Diff + IntDiff)^2;
      parameter Real r = 10.85/100; // meters
      parameter Real m = 6.8; // kg
      Real mu; // coefficient of friction, computed below
      constant Real g = 9.8; // gravity
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
      der(w_x)*cos(phi)+der(w_y)*sin(phi) = (1/I_x_p)*(m*r*(der(v_y)*cos(phi)+der(v_x)*sin(phi))+(I_y_p - I_z_p)*(-w_x*sin(phi)+w_y*cos(phi))*w_z);
//Equation 13
      -der(w_x)*sin(phi)+der(w_y)*cos(phi) = (1/I_y_p)*(m*r*(-der(v_y)*sin(phi)-der(v_x)*cos(phi))+(I_z_p - I_x_p)*(w_x*cos(phi)+w_y*sin(phi))*w_z);
//Equation 14
      der(w_z) = ((I_x_p - I_y_p)/I_z_p)*(w_x*cos(phi)+w_y*sin(phi))*(-w_x*sin(phi)+w_y*cos(phi));
//Equation 15
      der(phi) = der(atan(w_y/w_x));
      
      yFeet = y * 3.28084;
      xBoard = 39 - (x * 39.3701)/1.07; //Start numbering from right to left
  
      if y < 12.192 then
        mu = 0.04; // 0.04
      else
        mu = 0.2; // 0.2
      end if;
  annotation(
      experiment(StartTime = 0, StopTime = 3, Tolerance = 1e-06, Interval = 0.01));
  end Ball;
  
end Bowling;
