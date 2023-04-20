#include "StableSolver.h"
#include <cstdlib>
#include <array>
#include <cmath>

StableSolver::StableSolver(const Settings& aSetting)
{
    mySettings = aSetting;
    myTimeStep = 1.0f;

    const auto gridSize = GetGridSize();

    myVx.resize(gridSize);
    myVy.resize(gridSize);
    myVx0.resize(gridSize);
    myVy0.resize(gridSize);
    myD.resize(gridSize);
    myD0.resize(gridSize);
    myPx.resize(gridSize);
    myPy.resize(gridSize);
    myDiv.resize(gridSize);
    myP.resize(gridSize);
    myVort.resize(gridSize);
    myVortFx.resize(gridSize);
    myVortFy.resize(gridSize);

    for(int i=0; i<SIDE; i++)
    {
        for(int j=0; j<SIDE; j++)
        {
            myPx[IDX_1D(i, j)] = (float)i+0.5f;
            myPy[IDX_1D(i, j)] = (float)j+0.5f;
        }
    }
}

StableSolver::~StableSolver() {}

void StableSolver::ClearBuffer(const ClearMode aClearMode)
{
    switch (aClearMode) {
        case ClearMode::Internal: {
            std::fill(myD0.begin(), myD0.end(), 0);
            std::fill(myVx0.begin(), myVx0.end(), 0);
            std::fill(myVy0.begin(), myVy0.end(), 0);
            break;
        }
        case ClearMode::Sources: {
            std::fill(myD.begin(), myD.end(), 0);
            std::fill(myVx.begin(), myVx.end(), 0);
            std::fill(myVy.begin(), myVy.end(), 0);
            break;
        }
        default: break;
    }
}

void StableSolver::SetBoundary(int flag, std::vector<float>& anOutValue)
{
    static const std::array<std::pair<float, float>, 3> multipliers = {
            std::make_pair(1.f, 1.f),
            std::make_pair(1.f, -1.f),
            std::make_pair(-1.f, 1.f)
    };

    const auto & m = multipliers[flag];

    for(int i=1; i<SIDE-1; i++)
    {
        anOutValue[IDX_1D(i, 0)] = m.first * anOutValue[IDX_1D(i, 1)];
        anOutValue[IDX_1D(i, SIDE-1)] = m.first * anOutValue[IDX_1D(i, SIDE-2)];
        anOutValue[IDX_1D(0, i)] = m.second * anOutValue[IDX_1D(1, i)];
        anOutValue[IDX_1D(SIDE-1, i)] = m.second * anOutValue[IDX_1D(SIDE-2, i)];
    }

    anOutValue[IDX_1D(0, 0)] = (anOutValue[IDX_1D(0, 1)]+anOutValue[IDX_1D(1, 0)])/2;
    anOutValue[IDX_1D(SIDE-1, 0)] = (anOutValue[IDX_1D(SIDE-2, 0)]+anOutValue[IDX_1D(SIDE-1, 1)])/2;
    anOutValue[IDX_1D(0, SIDE-1)] = (anOutValue[IDX_1D(0, SIDE-2)]+anOutValue[IDX_1D(1, SIDE-1)])/2;
    anOutValue[IDX_1D(SIDE-1, SIDE-1)] = (anOutValue[IDX_1D(SIDE-2, SIDE-1)]+anOutValue[IDX_1D(SIDE-1, SIDE-2)])/2;
}

void StableSolver::Projection()
{
    for(int i=1; i<=SIDE-2; i++)
    {
        for(int j=1; j<=SIDE-2; j++)
        {
            myDiv[IDX_1D(i, j)] = 0.5f * (
                    myVx[IDX_1D(i+1, j)]-
                    myVx[IDX_1D(i-1, j)]+
                    myVy[IDX_1D(i, j+1)]-
                    myVy[IDX_1D(i, j-1)]);
            myP[IDX_1D(i, j)] = 0.0f;;
        }
    }
    SetBoundary(0, myDiv);
    SetBoundary(0, myP);

    constexpr auto projectionIterations = 20;
    constexpr auto factor = 1.f/4.f;
    for(int k=0; k<projectionIterations; k++)
    {
        for(int i=1; i<=SIDE-2; i++)
        {
            for(int j=1; j<=SIDE-2; j++)
            {
                myP[IDX_1D(i, j)] = (myP[IDX_1D(i+1, j)] +
                        myP[IDX_1D(i-1, j)] +
                        myP[IDX_1D(i, j+1)] +
                        myP[IDX_1D(i, j-1)] -
                        myDiv[IDX_1D(i, j)]) * factor;
            }
        }
        SetBoundary(0, myP);
    }

    for(int i=1; i<=SIDE-2; i++)
    {
        for(int j=1; j<=SIDE-2; j++)
        {
            myVx[IDX_1D(i, j)] -= 0.5f*(myP[IDX_1D(i+1, j)] - myP[IDX_1D(i-1, j)]);
            myVy[IDX_1D(i, j)] -= 0.5f*(myP[IDX_1D(i, j+1)] - myP[IDX_1D(i, j-1)]);
        }
    }

    SetBoundary(1, myVx);
    SetBoundary(2, myVy);
}

void StableSolver::Advection(std::vector<float>& value, const std::vector<float>& value0, const std::vector<float>& u, const std::vector<float>& v, int flag)
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
            oldX = myPx[IDX_1D(i, j)] - u[IDX_1D(i, j)]*myTimeStep;
            oldY = myPy[IDX_1D(i, j)] - v[IDX_1D(i, j)]*myTimeStep;
            
            if(oldX < 1.0f) oldX = 1.0f;
            if(oldX > SIDE-1) oldX = SIDE-1;
            if(oldY < 1.0f) oldY = 1.0f;
            if(oldY > SIDE-1) oldY = SIDE-1;

            i0 = (int)(oldX-0.5f);
            j0 = (int)(oldY-0.5f);
            i1 = i0+1;
            j1 = j0+1;
            
            wL = myPx[IDX_1D(i1, j0)]-oldX;
            wR = 1.0f-wL;
            wB = myPy[IDX_1D(i0, j1)]-oldY;
            wT = 1.0f-wB;

            value[IDX_1D(i, j)] = wB*(wL*value0[IDX_1D(i0, j0)]+wR*value0[IDX_1D(i1, j0)])+
                                wT*(wL*value0[IDX_1D(i0, j1)]+wR*value0[IDX_1D(i1, j1)]);
        }
    }
    
    SetBoundary(flag, value);
}

void StableSolver::Diffusion(std::vector<float>& value, const std::vector<float>& value0, float rate, int flag)
{
    const auto a = rate*myTimeStep;
    const auto factor = 1.f / (4.0f * a + 1.0f);

    constexpr auto diffusionIterations = 20;
    for(int k=0; k<diffusionIterations; k++)
    {
        for(int i=1; i<=SIDE-2; i++)
        {
            for(int j=1; j<=SIDE-2; j++)
            {
                value[IDX_1D(i, j)] = (value0[IDX_1D(i, j)] + a *
                        (value[IDX_1D(i+1, j)] +
                        value[IDX_1D(i-1, j)] +
                        value[IDX_1D(i, j+1)] +
                        value[IDX_1D(i, j-1)])) * factor;
            }
        }
        SetBoundary(flag, value);
    }
}

void StableSolver::VortexConfinement()
{
    for(int i=1; i<=SIDE-2; i++)
    {
        for(int j=1; j<=SIDE-2; j++)
        {
            myVort[IDX_1D(i, j)] = 0.5f * (myVy[IDX_1D(i+1, j)] - myVy[IDX_1D(i-1, j)] - myVx[IDX_1D(i, j+1)] + myVx[IDX_1D(i, j-1)]);
        }
    }

    SetBoundary(0, myVort);

    for(int i=1; i < SIDE-1; ++i)
    {
        for(int j=1; j < SIDE-1; ++j)
        {
            const auto currentGradVortX = 0.5f * (std::abs(myVort[IDX_1D(i+1, j)]) - std::abs(myVort[IDX_1D(i-1, j)]));
            const auto currentGradVortY = 0.5f * (std::abs(myVort[IDX_1D(i, j+1)]) - std::abs(myVort[IDX_1D(i, j-1)]));
            const auto currentLenGrad = std::sqrt(std::pow(currentGradVortX, 2) + std::pow(currentGradVortY, 2));

            constexpr auto zeroDistance = 0.01f;
            if(std::fabs(currentLenGrad) < zeroDistance)
            {
                myVortFx[IDX_1D(i, j)] = 0;
                myVortFy[IDX_1D(i, j)] = 0;
            }
            else
            {
                const auto factor = 1.f / currentLenGrad;
                myVortFx[IDX_1D(i, j)] = factor * currentGradVortX;
                myVortFy[IDX_1D(i, j)] = factor * currentGradVortY;
            }
        }
    }

    SetBoundary(0, myVortFx);
    SetBoundary(0, myVortFy);

    for(int i=1; i<=SIDE-2; i++)
    {
        for(int j=1; j<=SIDE-2; j++)
        {
            myVx[IDX_1D(i, j)] += mySettings.vorticity * (myVortFy[IDX_1D(i, j)] * myVort[IDX_1D(i, j)]);
            myVy[IDX_1D(i, j)] -= mySettings.vorticity * (myVortFx[IDX_1D(i, j)] * myVort[IDX_1D(i, j)]);
        }
    }

    SetBoundary(1, myVx);
    SetBoundary(2, myVy);
}

void StableSolver::AddSources()
{
    for(auto idx = 0; idx < SIDE*SIDE; ++idx)
    {
        myVx[idx] += myVx0[idx];
        myVy[idx] += myVy0[idx];
        myD[idx] += myD0[idx];
    }

    SetBoundary(0, myD);
    SetBoundary(1, myVx);
    SetBoundary(2, myVy);

}

void StableSolver::VelocityStep()
{
    if(mySettings.viscosity > 0.0f)
    {
        std::swap(myVx0, myVx);
        std::swap(myVy0, myVy);
        Diffusion(myVx, myVx0, mySettings.viscosity, 1);
        Diffusion(myVy, myVy0, mySettings.viscosity, 2);
    }

    Projection();

    std::swap(myVx0, myVx);
    std::swap(myVy0, myVy);
    Advection(myVx, myVx0, myVx0, myVy0, 1);
    Advection(myVy, myVy0, myVx0, myVy0, 2);

    Projection();
}

void StableSolver::DensityStep()
{
    if(mySettings.diffusion > 0.0f)
    {
        std::swap(myD0, myD);
        Diffusion(myD, myD0, mySettings.diffusion, 0);
    }

    std::swap(myD0, myD);
    Advection(myD, myD0, myVx, myVy, 0);
}
