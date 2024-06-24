#include <GL/glut.h>
#include <cmath>

// Function to draw a sphere
void drawSphere(double radius, int slices, int stacks) {
    glutSolidSphere(radius, slices, stacks);
}

// Function to draw a line based on Euler angles
void drawLine(double length, double yaw, double pitch, double roll) {
    // Convert Euler angles to radians
    yaw = yaw * M_PI / 180.0;
    pitch = pitch * M_PI / 180.0;
    roll = roll * M_PI / 180.0;

    // Calculate the endpoint of the line based on Euler angles
    double x = length * cos(yaw) * cos(pitch);
    double y = length * sin(pitch);
    double z = length * sin(yaw) * cos(pitch);

    // Draw the line
    glBegin(GL_LINES);
    glVertex3f(0.0, 0.0, 0.0); // Start point
    glVertex3f(x, y, z); // End point
    glEnd();
}

// Function to display OpenGL scene
void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    // Set camera position
    gluLookAt(5.0, 5.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

    // Draw the sphere
    glColor3f(1.0, 1.0, 1.0); // Set sphere color to white
    drawSphere(1.0, 50, 50);

    // Draw a line based on Euler angles
    glColor3f(1.0, 0.0, 0.0); // Set line color to red
    drawLine(1.5, 45.0, 30.0, 0.0); // Length and Euler angles (yaw, pitch, roll)

    glutSwapBuffers();
}

// Main function
int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(800, 600);
    glutCreateWindow("Sphere and Line");

    glEnable(GL_DEPTH_TEST);

    glutDisplayFunc(display);

    glutMainLoop();

    return 0;
}
