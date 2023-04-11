#include "StableSolver.h"
#include <cstdlib>
#include <cstring>
#include <array>
#include <cmath>

StableSolver::StableSolver(const Settings& aSetting)
{
    mySettings = aSetting;
    timeStep = 1.0f;

    const auto gridSize = GetGridSize();

    vx = (float *)malloc(sizeof(float)*gridSize);
    vy = (float *)malloc(sizeof(float)*gridSize);
    vx0 = (float *)malloc(sizeof(float)*gridSize);
    vy0 = (float *)malloc(sizeof(float)*gridSize);
    d = (float *)malloc(sizeof(float)*gridSize);
    d0 = (float *)malloc(sizeof(float)*gridSize);
    px = (float *)malloc(sizeof(float)*gridSize);
    py = (float *)malloc(sizeof(float)*gridSize);
    div = (float *)malloc(sizeof(float)*gridSize);
    p = (float *)malloc(sizeof(float)*gridSize);

    //vorticity confinement
    vort = (float *)malloc(sizeof(float)*gridSize);
    absVort = (float *)malloc(sizeof(float)*gridSize);
    gradVortX = (float *)malloc(sizeof(float)*gridSize);
    gradVortY = (float *)malloc(sizeof(float)*gridSize);
    lenGrad = (float *)malloc(sizeof(float)*gridSize);
    vcfx = (float *)malloc(sizeof(float)*gridSize);
    vcfy = (float *)malloc(sizeof(float)*gridSize);

    for(int i=0; i<SIDE; i++)
    {
        for(int j=0; j<SIDE; j++)
        {
            px[IDX_1D(i, j)] = (float)i+0.5f;
            py[IDX_1D(i, j)] = (float)j+0.5f;
        }
    }
}

StableSolver::~StableSolver()
{
    free(vx);
    free(vy);
    free(vx0);
    free(vy0);
    free(d);
    free(d0);
    free(px);
    free(py);
    free(div);
    free(p);

    //vorticity confinement
    free(vort);
    free(absVort);
    free(gradVortX);
    free(gradVortY);
    free(lenGrad);
    free(vcfx);
    free(vcfy);
}

void StableSolver::ClearBuffer(const ClearMode aClearMode)
{
    const auto gridSize = GetGridSize();

    switch (aClearMode) {
        case ClearMode::Internal: {
            memset(d0, 0, sizeof(float) * gridSize);
            memset(vx0, 0, sizeof(float) * gridSize);
            memset(vy0, 0, sizeof(float) * gridSize);
            break;
        }
        case ClearMode::Sources: {
            memset(d, 0, sizeof(float)*gridSize);
            memset(vx, 0, sizeof(float)*gridSize);
            memset(vy, 0, sizeof(float)*gridSize);
            break;
        }
        default: break;
    }
}

void StableSolver::SetBoundary(float *value, int flag)
{
    static const std::array<std::pair<float, float>, 3> multipliers = {
            std::make_pair(1.f, 1.f),
            std::make_pair(1.f, -1.f),
            std::make_pair(-1.f, 1.f)
    };

    const auto & m = multipliers[flag];

    for(int i=1; i<SIDE-1; i++)
    {
        value[IDX_1D(i, 0)] = m.first * value[IDX_1D(i, 1)];
        value[IDX_1D(i, SIDE-1)] = m.first * value[IDX_1D(i, SIDE-2)];
        value[IDX_1D(0, i)] = m.second * value[IDX_1D(1, i)];
        value[IDX_1D(SIDE-1, i)] = m.second * value[IDX_1D(SIDE-2, i)];
    }

    value[IDX_1D(0, 0)] = (value[IDX_1D(0, 1)]+value[IDX_1D(1, 0)])/2;
    value[IDX_1D(SIDE-1, 0)] = (value[IDX_1D(SIDE-2, 0)]+value[IDX_1D(SIDE-1, 1)])/2;
    value[IDX_1D(0, SIDE-1)] = (value[IDX_1D(0, SIDE-2)]+value[IDX_1D(1, SIDE-1)])/2;
    value[IDX_1D(SIDE-1, SIDE-1)] = (value[IDX_1D(SIDE-2, SIDE-1)]+value[IDX_1D(SIDE-1, SIDE-2)])/2;
}

void StableSolver::Projection()
{
    for(int i=1; i<=SIDE-2; i++)
    {
        for(int j=1; j<=SIDE-2; j++)
        {
            div[IDX_1D(i, j)] = 0.5f * (
                    vx[IDX_1D(i+1, j)]-
                    vx[IDX_1D(i-1, j)]+
                    vy[IDX_1D(i, j+1)]-
                    vy[IDX_1D(i, j-1)]);
            p[IDX_1D(i, j)] = 0.0f;;
        }
    }
    SetBoundary(div, 0);
    SetBoundary(p, 0);

    constexpr auto projectionIterations = 20;
    for(int k=0; k<projectionIterations; k++)
    {
        for(int i=1; i<=SIDE-2; i++)
        {
            for(int j=1; j<=SIDE-2; j++)
            {
                p[IDX_1D(i, j)] = (p[IDX_1D(i+1, j)]+p[IDX_1D(i-1, j)]+p[IDX_1D(i, j+1)]+p[IDX_1D(i, j-1)]-div[IDX_1D(i, j)])/4.0f;
            }
        }
        SetBoundary(p, 0);
    }

    for(int i=1; i<=SIDE-2; i++)
    {
        for(int j=1; j<=SIDE-2; j++)
        {
            vx[IDX_1D(i, j)] -= 0.5f*(p[IDX_1D(i+1, j)]-p[IDX_1D(i-1, j)]);
            vy[IDX_1D(i, j)] -= 0.5f*(p[IDX_1D(i, j+1)]-p[IDX_1D(i, j-1)]);
        }
    }
    SetBoundary(vx, 1);
    SetBoundary(vy, 2);
}

void StableSolver::Advection(float *value, const float *value0, const float *u, const float *v, int flag)
{
    float oldX;
    float oldY;
    int i0;
    int i1;
    int j0;
    int j1;
    float wL;
    float wR;
    float wB;
    float wT;

    for(int i=1; i<=SIDE-2; i++)
    {
        for(int j=1; j<=SIDE-2; j++)
        {
            oldX = px[IDX_1D(i, j)] - u[IDX_1D(i, j)]*timeStep;
            oldY = py[IDX_1D(i, j)] - v[IDX_1D(i, j)]*timeStep;
            
            if(oldX < 1.0f) oldX = 1.0f;
            if(oldX > SIDE-1) oldX = SIDE-1;
            if(oldY < 1.0f) oldY = 1.0f;
            if(oldY > SIDE-1) oldY = SIDE-1;

            i0 = (int)(oldX-0.5f);
            j0 = (int)(oldY-0.5f);
            i1 = i0+1;
            j1 = j0+1;
            
            wL = px[IDX_1D(i1, j0)]-oldX;
            wR = 1.0f-wL;
            wB = py[IDX_1D(i0, j1)]-oldY;
            wT = 1.0f-wB;

            value[IDX_1D(i, j)] = wB*(wL*value0[IDX_1D(i0, j0)]+wR*value0[IDX_1D(i1, j0)])+
                                wT*(wL*value0[IDX_1D(i0, j1)]+wR*value0[IDX_1D(i1, j1)]);
        }
    }
    
    SetBoundary(value, flag);
}

void StableSolver::Diffusion(float *value, const float *value0, float rate, int flag)
{
    float a = rate*timeStep;

    constexpr auto diffusionIterations = 20;
    for(int k=0; k<diffusionIterations; k++)
    {
        for(int i=1; i<=SIDE-2; i++)
        {
            for(int j=1; j<=SIDE-2; j++)
            {
                value[IDX_1D(i, j)] = (value0[IDX_1D(i, j)]+a*(value[IDX_1D(i+1, j)]+value[IDX_1D(i-1, j)]+value[IDX_1D(i, j+1)]+value[IDX_1D(i, j-1)])) / (4.0f*a+1.0f);
            }
        }
        SetBoundary(value, flag);
    }
}

void StableSolver::VortexConfinement()
{
    for(int i=1; i<=SIDE-2; i++)
    {
        for(int j=1; j<=SIDE-2; j++)
        {
            vort[IDX_1D(i, j)] = 0.5f*(vy[IDX_1D(i+1, j)]-vy[IDX_1D(i-1, j)]-vx[IDX_1D(i, j+1)]+vx[IDX_1D(i, j-1)]);
            if(vort[IDX_1D(i, j)] >= 0.0f) absVort[IDX_1D(i, j)] = vort[IDX_1D(i, j)];
            else absVort[IDX_1D(i, j)] = -vort[IDX_1D(i, j)];
        }
    }
    SetBoundary(vort, 0);
    SetBoundary(absVort, 0);

    for(int i=1; i<=SIDE-2; i++)
    {
        for(int j=1; j<=SIDE-2; j++)
        {
            gradVortX[IDX_1D(i, j)] = 0.5f*(absVort[IDX_1D(i+1, j)]-absVort[IDX_1D(i-1, j)]);
            gradVortY[IDX_1D(i, j)] = 0.5f*(absVort[IDX_1D(i, j+1)]-absVort[IDX_1D(i, j-1)]);
            lenGrad[IDX_1D(i, j)] = sqrt(gradVortX[IDX_1D(i, j)]*gradVortX[IDX_1D(i, j)]+gradVortY[IDX_1D(i, j)]*gradVortY[IDX_1D(i, j)]);
            if(lenGrad[IDX_1D(i, j)] < 0.01f)
            {
                vcfx[IDX_1D(i, j)] = 0.0f;
                vcfy[IDX_1D(i, j)] = 0.0f;
            }
            else
            {
                vcfx[IDX_1D(i, j)] = gradVortX[IDX_1D(i, j)] / lenGrad[IDX_1D(i, j)];
                vcfy[IDX_1D(i, j)] = gradVortY[IDX_1D(i, j)] / lenGrad[IDX_1D(i, j)];
            }
        }
    }
    SetBoundary(vcfx, 0);
    SetBoundary(vcfy, 0);

    for(int i=1; i<=SIDE-2; i++)
    {
        for(int j=1; j<=SIDE-2; j++)
        {
            vx[IDX_1D(i, j)] += mySettings.vorticity * (vcfy[IDX_1D(i, j)] * vort[IDX_1D(i, j)]);
            vy[IDX_1D(i, j)] -= mySettings.vorticity * (vcfx[IDX_1D(i, j)] * vort[IDX_1D(i, j)]);
        }
    }

    SetBoundary(vx, 1);
    SetBoundary(vy, 2);
}

void StableSolver::AddSources()
{
    for(auto idx = 0; idx < SIDE*SIDE; ++idx)
    {
        vx[idx] += vx0[idx];
        vy[idx] += vy0[idx];
        d[idx] += d0[idx];
    }

    SetBoundary(d, 0);
    SetBoundary(vx, 1);
    SetBoundary(vy, 2);

}

void StableSolver::VelocityStep()
{
    if(mySettings.viscosity > 0.0f)
    {
        std::swap(vx0, vx);
        std::swap(vy0, vy);
        Diffusion(vx, vx0, mySettings.viscosity, 1);
        Diffusion(vy, vy0, mySettings.viscosity, 2);
    }

    Projection();

    std::swap(vx0, vx);
    std::swap(vy0, vy);
    Advection(vx, vx0, vx0, vy0, 1);
    Advection(vy, vy0, vx0, vy0, 2);

    Projection();
}

void StableSolver::DensityStep()
{
    if(mySettings.diffusion > 0.0f)
    {
        std::swap(d0, d);
        Diffusion(d, d0, mySettings.diffusion, 0);
    }

    std::swap(d0, d);
    Advection(d, d0, vx, vy, 0);
}
