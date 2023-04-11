#include "Visualizer.h"
#include "StableSolver.h"
#include <GLUT/glut.h>

void VelocityVisualizer::operator()()
{
    glColor3f(1.f, 0.f, 0.f);
    glLineWidth(1.0f);

    glBegin(GL_LINES);
    for(int i=0; i<mySolver->GetGridSize(); i++)
    {
        const auto pressure = mySolver->GetPressureAt(i);
        const auto velocity = mySolver->GetVelocityAt(i);
        glVertex2f(pressure.first, pressure.second);
        glVertex2f(pressure.first+velocity.first*10.0f, pressure.second+velocity.second*10.0f);
    }
    glEnd ();
}

void DensityVisualizer::operator()()
{
    float x, y, d00, d01, d10, d11;

    glBegin(GL_QUADS);
    for(int i=1; i<mySolver->GetSide()-1; i++)
    {
        x = (float)i;
        for(int j=1; j<mySolver->GetSide()-1; j++)
        {
            y = (float)j;

            d00 = mySolver->GetDensityAt(i, j);
            d01 = mySolver->GetDensityAt(i, j+1);
            d10 = mySolver->GetDensityAt(i+1, j);
            d11 = mySolver->GetDensityAt(i+1, j+1);

            glColor3f(1.0f, 1.0f-d00, 1.0f-d00); glVertex2f(x, y);
            glColor3f(1.0f, 1.0f-d10, 1.0f-d10); glVertex2f(x+1.0f, y);
            glColor3f(1.0f, 1.0f-d11, 1.0f-d11); glVertex2f(x+1.0f, y+1.0f);
            glColor3f(1.0f, 1.0f-d01, 1.0f-d01); glVertex2f(x, y+1.0f);
        }
    }
    glEnd();
}
