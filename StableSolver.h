#ifndef GAMEPHYSICSPROJECT_STABLESOLVER_H
#define GAMEPHYSICSPROJECT_STABLESOLVER_H

#include <vector>
#include <cassert>

#define SIDE 130
#define IDX_1D(x, y) ((y)*(SIDE)+(x))

class StableSolver
{
public:

    enum class State
    {
        Stopped,
        Running
    };

    struct Settings
    {
        State myState {State::Stopped};
        float viscosity {0};
        float diffusion {0};
        float vorticity {0};
    };

    enum class ClearMode
    {
        Sources,
        Internal
    };

    StableSolver(const Settings& aSetting);
    ~StableSolver();

    void ClearBuffer(const ClearMode aClearMode = ClearMode::Internal);

    void ToggleExecution()
    {
        mySettings.myState = static_cast<State>(1 ^ static_cast<int>(mySettings.myState));
    }

    void UpdateSettings(const Settings& aSettings)
    {
        mySettings = aSettings;
    }

    const Settings& GetSettings() const
    {
        return mySettings;
    }

    void SimulationStep()
    {
        if(mySettings.myState == State::Running)
        {
            AddSources();
            VortexConfinement();
            VelocityStep();
            DensityStep();
        }
    }

    int GetSide() const { return SIDE; }
    int GetGridSize() const { return std::pow(SIDE, 2); }

    std::pair<float, float> GetVelocityAt(int index) const
    {
        assert(index > -1 && index < GetGridSize());
        return std::make_pair(myVx[index], myVy[index]);
    }

    std::pair<float, float> GetPressureAt(int index) const
    {
        assert(index > -1 && index < GetGridSize());
        return std::make_pair(myPx[index], myPy[index]);
    }

    float GetDensityAt(int i, int j) const
    {
        constexpr auto factor = 1.f / 4.f;

        assert(i > -1 && i < SIDE && j > -1 && j < SIDE);

        return (myD[IDX_1D(i-1, j-1)] +
                myD[IDX_1D(i, j-1)] +
                myD[IDX_1D(i-1, j)] +
                myD[IDX_1D(i, j)]) * factor;
    }

    void SetVelocitySource(int i, int j, const std::pair<float, float>& aValue) const
    {
        assert(i > -1 && i < SIDE && j > -1 && j < SIDE);

        if(i == 0 || i == SIDE-1 || j == 0 || j == SIDE-1)
        {
            return;
        }

        const auto index = IDX_1D(i, j);
        myVx0[index] = aValue.first;
        myVy0[index] = aValue.second;
    }

    void SetDensitySource(int i, int j, const float aValue) const
    {
        assert(i > -1 && i < SIDE && j > -1 && j < SIDE);

        if(i == 0 || i == SIDE-1 || j == 0 || j == SIDE-1)
        {
            return;
        }

        myD0[IDX_1D(i, j)] = aValue;
    }


private:
    void SetBoundary(float *value, int flag);
    void Projection();
    void Advection(float *value, const float *value0, const float *u, const float *v, int flag);
    void Diffusion(float *value, const float *value0, float rate, int flag);
    void VortexConfinement();
    void AddSources();
    void VelocityStep();
    void DensityStep();


    Settings mySettings;
    float myTimeStep;

    float* myVx;
    float* myVy;
    float* myVx0;
    float* myVy0;
    float * myD;
    float * myD0;
    float * myPx;
    float * myPy;
    float * myDiv;
    float * myP;

    float * myVort;
    float * myVortFx;
    float * myVortFy;
};

#endif
