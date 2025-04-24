encapsulated package Bowling
  type Velocity = Real(unit = "m/s");
  type AngularVelocity = Real(unit = "rad/s");
  type Inertia = Real(unit = "I");
  type Force = Real(unit="N");
  model Bowl
    //
  
  
  end Bowl;
  model Ball

  

  
  equation
  v_b_x = v_x - (r * w_y);
  v_b_y = v_y + (r * w_x);
  
  end Ball;
  
  //Rotation of Ball
  package Rotational
    type Torque = Real(unit = "N.m");
    type Angle = Real(unit = "rad");
    type AngularVelocity = Real(unit = "rad/s");
    
    connector Contact
      Angle th;
      flow Torque tau;
      AngularVelocity w_x = -30;
      AngularVelocity w_y = -30;
      AngularVelocity w_z = 10;
    end Contact;
    
    model Ball
      Rotational.Inertia inertia;
      Rotational.friction friction;

    equation
      connect(friction.contact,inertia.contact);
    end Ball;
    

  
    model Friction
      AngularVelocity w;
      parameter Real b = 1;
      Rotational.Contact contact;
    equation
      contact.tau = b * w;
      der(contact.th) = w;
    end Friction;
    
  end Rotational;  
  
  //Movement on Lane
  package Translational
    type Position = Real(unit="m");
    type Distance = Real(unit="m");
    type Velocity = Real(unit="m/s");
    type Force = Real(unit="N");

    connector Contact
      Position x;
      Position y;
      Velocity vx;
      Velocity vy;
      flow Force f;
    end Contact;
    
    model Ball
      Translational.Contact contact;
      parameter Real r = 10.85;
      parameter Real m = 6.8;
    equation
      x = contact.x;
      y = contact.y;
      v = der(x);
      m*der(v) = contact.f;
    end Ball;
      
  end Translational;
  
  package RotoTranslational
  
    model Ball
      parameter Velocity v_x = 8;
      parameter Velocity v_y = 8;
      parameter AngularVelocity w_x = -30;
      parameter AngularVelocity w_y = -30;
      parameter AngularVelocity w_z = 10;
      Translational.Ball translational;
      Rotational.Ball rotational;
      Translational.Contact tContact;
      Rotational.Contact rContact;
    equation

    end Ball;
      
    model Inertia
      AngularVelocity w;
      parameter Real RG = 6.35;
      parameter Real Diff = 0.1;
      parameter Real IntDiff = 0;
      parameter Inertia I_y = m * (RG)^2;
      parameter Inertia I_x = m * (RG + Diff)^2;
      parameter Inertia I_z = m * (RG + Diff + IntDiff)^2;
      parameter Real r = 10.85;
      parameter Real m = 6.8;
      Rotational.Contact contact;
    equation
      //Equation 12
      der(contact.w_x)*cos(contact.th)+der(contact.w_y)*sin(contact.th) = 
      (1/I_x)*
      (m*r*(der(vy)*cos(contact.th)+der(vx)*sin(contact.th))
      +(I_y - I_z)*(-w_x*sin(contact.th)+contact.w_y*cos(contact.th)*contact.w_z));
      //Equation 13
      -der(contact.w_x)*sin(contact.th)+der(contact.w_y)*cos(contact.th) = 
      (1/I_y)*
      (m*r*(-der(vy)*sin(contact.th)-der(vx)*cos(contact.th))
      +(I_z - I_x)*(w_x*cos(contact.th)+w_y*sin(contact.th)*w_z));
      //Equation 14
      der(w_z) = ((I_x - I_y)/I_z)*
      (w_x*cos(contact.th)+w_y*sin(contact.th))*
      (-w_x*sin(contact.th)+w_y*cos(contact.th));
    end Inertia;
  
  end RotoTranslational;
end Bowling;
