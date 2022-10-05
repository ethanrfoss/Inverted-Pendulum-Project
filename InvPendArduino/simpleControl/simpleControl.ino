

#include <Thread.h>
#include <ThreadController.h>
#include <RobojaxBTS7960.h>

//Threads
Thread inLoop = Thread();
Thread outLoop = Thread();
ThreadController threads = ThreadController();

//Limit Switch Pins
#define limLPin 53
#define limRPin 52

//DC Motor Pins
#define RPWM 3 // define pin 3 for RPWM pin (output)
#define R_EN 4 // define pin 2 for R_EN pin (input)
#define R_IS 5 // define pin 5 for R_IS pin (output)

#define LPWM 6 // define pin 6 for LPWM pin (output)
#define L_EN 7 // define pin 7 for L_EN pin (input)
#define L_IS 8 // define pin 8 for L_IS pin (output)
#define CW 1 //do not change
#define CCW 0 //do not change

//DC Motor Object
RobojaxBTS7960 motor(R_EN, RPWM, R_IS, L_EN, LPWM, L_IS, 1);

//Encoder Pins
#define w1 18
#define g1 19
#define w2 20
#define g2 21

//Encoder Signals
long counter1 = 0;
long counter2 = 0;

//Loop Signals
float ex;
float et;
float tCom;
float uCom;
float lastTime;

//System References
float centerPos;
float pendulumDownAngle;
float pendulumUpAngle;

//System Parameters
float r = .021;

void setup() {
  Serial.begin(115200);

  //Motor
  motor.begin();
  
  // Threads
  inLoop.onRun(innerLoop);
  inLoop.setInterval(20);
  outLoop.onRun(outerLoop);
  outLoop.setInterval(100);
  threads.add(&inLoop);
  threads.add(&outLoop);

  // Limit Switch PinsSetip
  pinMode(limLPin,INPUT_PULLUP);
  pinMode(limRPin,INPUT_PULLUP);

  // Encoder Pin Setup
  pinMode(w1, INPUT_PULLUP);
  pinMode(g1, INPUT_PULLUP);
  attachInterrupt(digitalPinToInterrupt(g1), encoder1PinChange, CHANGE);
  pinMode(w2, INPUT_PULLUP);
  pinMode(g2, INPUT_PULLUP);
  attachInterrupt(digitalPinToInterrupt(g2), encoder2PinChange, CHANGE);
  
  getDownAngle();
  centerGantry();
  lastTime = micros();
}

void loop() {
  Serial.println(getPendAng()-pendulumUpAngle,4);
  if(digitalRead(limRPin) == 0 || digitalRead(limLPin) == 0){
    motor.stop();
  }
  else if(getPendAng() > pendulumUpAngle + .1 || getPendAng() < pendulumUpAngle - .1){
    motor.stop();
  }
  else{
    threads.run();
  }
}

void getDownAngle(){

    float minAngle = getPendAng();
    float maxAngle = getPendAng();

    float startTime = micros();
    float ang;

    while(micros() <= startTime + 5000000){
      ang = getPendAng();
      if(minAngle>ang){
        minAngle = ang;
      }
      if(maxAngle<ang){
        maxAngle = ang;
      }
    }

    pendulumDownAngle = (maxAngle+minAngle)/2;
    pendulumUpAngle = pendulumDownAngle + 3.14159;
}

void centerGantry(){

  while(digitalRead(limRPin) == 1){
    motor.rotate(15,CW);
  }
  motor.stop();
  float rightPos = getPos();

  while(digitalRead(limLPin) == 1){
    motor.rotate(15,CCW);
  }
  motor.stop();
  float leftPos = getPos();

  centerPos = (rightPos+leftPos)/2;

  float startTime = micros();
  float integral = 0;
  float oldTime = micros();
  float newTime = micros();

  Serial.println(getPos()-centerPos);
  Serial.println(1.25*(getPos()-centerPos));
  
  while(newTime-startTime<=3000000){
    oldTime = newTime;
    newTime = micros();
    integral = integral + (centerPos-getPos())*(newTime-oldTime)/(1000000);
    motor.rotate((1.25*(centerPos-getPos())+.25*integral)*100,CW);
  }
  motor.stop();
}

void innerLoop(){

  //float etNew = tCom-(getPendAng()-pendulumUpAngle);
  //float uComNew = .9985*uCom + 7.846*etNew -7.422*et;

  //et = etNew;
  //uCom = uComNew;

  float etNew = tCom-(getPendAng()-pendulumUpAngle);
  float uComNew = .997*uCom + 8.058*etNew - 7.21*et;

  et = etNew;
  uCom = uComNew;

  if(uCom >= 0){
    motor.rotate(min(uCom,1)*100,CW);
  }
  else{
    motor.rotate(min(-uCom,1)*100,CCW);
  }

  Serial.print("Input Command: ");
  Serial.print(uCom);
  Serial.print(" Theta Command: ");
  Serial.print(tCom);
  Serial.print(" Current Angle: ");
  Serial.print(getPendAng()-pendulumUpAngle);
  Serial.print(" Current Position: ");
  Serial.print(getPos());
  Serial.print(" Time Int: ");
  Serial.println((micros()-lastTime)/1000000);
  lastTime = micros();
  
}

void outerLoop(){

//  float exNew = -(centerPos - getPos());
//  float tComNew = .6065*tCom - .2424*exNew + .2354*ex;
//
//  ex = exNew;
//  tCom = tComNew;

  float exNew = (centerPos - getPos());
  float tComNew = .3679*tCom - .1976*exNew + .1863*ex;

  ex = exNew;
  tCom = tComNew;
  
}

float getPendAng(){
  float angle;
  angle = -counter1/1200.0*2*3.14159;
  return angle;
}

float getPos(){
  float pos;
  pos = -counter2/1200.0*2*3.14159*r;
  return pos;
}

void encoder1PinChange() {
counter1 += digitalRead(g1) == digitalRead(w1) ? -1 : 1;
}

void encoder2PinChange() {
counter2 += digitalRead(g2) != digitalRead(w2) ? -1 : 1;
}
