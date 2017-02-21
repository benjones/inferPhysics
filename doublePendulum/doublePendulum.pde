
float g = 10;
float dt = .01;
float m1 = 2;
float m2 = .05;
float l1 = 1;
float l2 = .75;

float theta1, theta2, theta1Dot, theta2Dot;

void setup() {
  size(600, 600);
  frameRate(60);
  background(255);

  theta1 = 3*PI/4;
  theta2 = .9*PI;
  theta1Dot = 0;
  theta2Dot = 0;
}

void drawPoint(float x, float y) {

  point(x*width*.25 + width*.5, -y*height*.25 + height*.5);
}

void drawLine(float x1, float y1, float x2, float y2) {
  line(x1*width*.25 + width*.5, -y1*height*.25 + height*.5, 
    x2*width*.25 + width*.5, -y2*height*.25 + height*.5);
}

void draw() {
  background(255);

  float a = (m1 + m2)*l1;
  float b = m2*l2*cos(theta1 - theta2);

  float c = -m2*l2*theta2Dot*theta2Dot*sin(theta1 - theta2) - 
    g*(m1 + m2)*sin(theta1);

  float d = m2*l1*cos(theta1 - theta2);
  float e = m2*l2;
  float f = m2*l1*theta1Dot*theta1Dot*sin(theta1 - theta2) - 
    m2*g*sin(theta2);

  float theta2DotDot = (f - c*d/a)/(e - b*d/a);
  float theta1DotDot = (c - b*theta2DotDot)/a;

  theta1Dot += dt*theta1DotDot;
  theta2Dot += dt*theta2DotDot;

  theta1Dot *= .999999;
  theta2Dot *= .999999;
  
  theta1 += dt*theta1Dot;
  theta2 += dt*theta2Dot;

  float x1 = l1*sin(theta1);
  float y1 = -l1*cos(theta1);
  float x2 = x1 + l2*sin(theta2);
  float y2 = y1 - l2*cos(theta2);

  stroke(0);
  strokeWeight(1);
  drawLine(0, 0, x1, y1);
  drawLine(x1, y1, x2, y2);
  strokeWeight(10);
  drawPoint(0, 0);
  stroke(255, 0, 0);
  drawPoint(x1, y1);
  stroke(0, 0, 255);
  drawPoint(x2, y2);
}