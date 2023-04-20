#include <GLUT/glut.h>
#include <cstdlib>
#include "Visualizer.h"
#include "StableSolver.h"
#include "externals/imgui/imgui.h"
#include "wrapper/imgui_impl_glut.h"
#include "wrapper/imgui_impl_opengl2.h"

#define WIDTH 800
#define HEIGHT 800

constexpr StableSolver::Settings defaultSettings
{
    StableSolver::State::Running,
    0.003f,
    0.01f,
    0.01f
};

StableSolver* solver = new StableSolver(defaultSettings);

int mouseDownEvents[] = {0, 0, 0};
int origin_x, origin_y, dest_x, dest_y;

void GetMouseInput()
{
    solver->ClearBuffer();

    if(mouseDownEvents[0] || mouseDownEvents[2])
    {
        const int size = solver->GetSide();
        const int x = (int)((float)(origin_x)/WIDTH*size);
        const int y = (int)((float)(HEIGHT - origin_y)/HEIGHT*size);

        if(x < 0 || x >= size || y < 0 || y >= size)
        {
            return;
        }

        if(mouseDownEvents[0])
        {
            solver->SetVelocitySource(x, y, std::make_pair(dest_x - origin_x, origin_y - dest_y));
        }
        else if(mouseDownEvents[2])
        {

            solver->SetDensitySource(x, y, 10.0f);
        }

        origin_x = dest_x;
        origin_y = dest_y;
    }
}

void OnKeyPressed(unsigned char key, int x, int y)
{
    switch(key)
    {
        case 27: // escape
            exit(0);
            break;
    }
}

void OnMousePressed(int button, int state, int x, int y)
{
    ImGuiIO& io = ImGui::GetIO();

    if(io.WantCaptureMouse)
    {
        io.AddMousePosEvent((float)x, (float)y);

        if(state == GLUT_DOWN || state == GLUT_UP)
        {
            switch (button) {
                case GLUT_LEFT_BUTTON:
                    io.AddMouseButtonEvent(0, state == GLUT_DOWN);
                    break;
                case GLUT_RIGHT_BUTTON:
                    io.AddMouseButtonEvent(1, state == GLUT_DOWN);
                    break;
                default:
                    break;
            }

            if(state == GLUT_UP)
            {
                mouseDownEvents[button] = state == GLUT_DOWN;
            }
        }
    }
    else
    {
        origin_x = x;
        origin_y = y;

        dest_x = x;
        dest_y = y;

        mouseDownEvents[button] = state == GLUT_DOWN;
    }
}

void OnMouseMoved(int x, int y)
{
    ImGuiIO& io = ImGui::GetIO();
    io.AddMousePosEvent((float)x, (float)y);

    dest_x = x;
    dest_y = y;
}

void ReshapeWindow (int, int)
{
    glutReshapeWindow (WIDTH, HEIGHT);
    ImGuiIO& io = ImGui::GetIO();
    io.DisplaySize = ImVec2(WIDTH, HEIGHT);
}

 void DisplayWindow()
 {
    static int visualizationTypeIndex = 1;

     const std::vector<Visualizer*> visualizers = {
             new DensityVisualizer{solver},
             new VelocityVisualizer{solver}
     };

     GetMouseInput();
    solver->SimulationStep();

     // Start the Dear ImGui frame
     ImGui_ImplOpenGL2_NewFrame();
     ImGui_ImplGLUT_NewFrame();

     ImGui::SetNextWindowSize(ImVec2(300, 240));
     ImGui::Begin("Settings");


     if (ImGui::CollapsingHeader("StableSolver"))
     {
         auto settings = solver->GetSettings();

         if(    ImGui::SliderFloat("Diffusion", &settings.diffusion, 0.0f, 1.0f) ||
                ImGui::SliderFloat("Viscosity", &settings.viscosity, 0.0f, 1.0f) ||
                ImGui::SliderFloat("Vorticity", &settings.vorticity, 0.0f, 1.0f))
         {
             solver->UpdateSettings(settings);
         }

         if(ImGui::Button("Reset"))
         {
             solver->ClearBuffer(StableSolver::ClearMode::Sources);
             solver->ClearBuffer(StableSolver::ClearMode::Internal);
         }

         if(ImGui::Button("Pause/Unpause"))
         {
             solver->ToggleExecution();
         }
     }

     if(ImGui::CollapsingHeader("Visualization"))
     {
         ImGui::RadioButton("Density", &visualizationTypeIndex, 0); ImGui::SameLine();
         ImGui::RadioButton("Velocity", &visualizationTypeIndex, 1);
     }

     ImGuiIO& io = ImGui::GetIO(); (void)io;
     ImGui::Text("AVG %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
     ImGui::End();

     // Rendering
     glViewport(0, 0, WIDTH, HEIGHT);
     glMatrixMode(GL_PROJECTION);
     glLoadIdentity ();
     gluOrtho2D(0.0f, (float)(solver->GetSide()), 0.0f, (float)(solver->GetSide()));
     glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
     glClear(GL_COLOR_BUFFER_BIT);

     visualizers[visualizationTypeIndex]->operator()();

     glColor3f(0.0f, 0.0f, 0.0f);
     glPointSize(1.0f);
     glBegin(GL_POINTS);
     for(int i=0; i<solver->GetGridSize(); i++)
     {
         const auto pressure = solver->GetPressureAt(i);
         glVertex2f(pressure.first, pressure.second);
     }
     glEnd();

     ImGui::Render();
     ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());

     glutSwapBuffers();
     glutPostRedisplay();
}

int main(int argc, char** argv)
{
    solver->ClearBuffer(StableSolver::ClearMode::Internal);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(WIDTH, HEIGHT);
    glutCreateWindow("2D Stable Fluid");

    glutKeyboardFunc(OnKeyPressed);
    glutMouseFunc(OnMousePressed);
    glutMotionFunc(OnMouseMoved);
    glutPassiveMotionFunc(OnMouseMoved);
    glutReshapeFunc(ReshapeWindow);
    glutIdleFunc(glutPostRedisplay);
    glutDisplayFunc(DisplayWindow);

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();

    // Setup Platform/Renderer backends
    ImGui_ImplGLUT_Init();
    ImGui_ImplOpenGL2_Init();

    ImGui_ImplGLUT_InstallFuncs();

    glutMainLoop();

    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGLUT_Shutdown();
    ImGui::DestroyContext();

    return 0;
}
