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

    void SetVelocitySource(int i, int j, const std::pair<float, float>& aValue)
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

    void SetDensitySource(int i, int j, const float aValue)
    {
        assert(i > -1 && i < SIDE && j > -1 && j < SIDE);

        if(i == 0 || i == SIDE-1 || j == 0 || j == SIDE-1)
        {
            return;
        }

        myD0[IDX_1D(i, j)] = aValue;
    }


private:
    void SetBoundary(int flag, std::vector<float>& anOutValue);
    void Projection();
    void Advection(std::vector<float>& value, const std::vector<float>& value0, const std::vector<float>& u, const std::vector<float>& v, int flag);
    void Diffusion(std::vector<float>& value, const std::vector<float>& value0, float rate, int flag);
    void VortexConfinement();
    void AddSources();
    void VelocityStep();
    void DensityStep();


    Settings mySettings;
    float myTimeStep;

    std::vector<float> myVx;
    std::vector<float> myVy;
    std::vector<float> myVx0;
    std::vector<float> myVy0;

    std::vector<float> myD;
    std::vector<float> myD0;
    std::vector<float> myPx;
    std::vector<float> myPy;
    std::vector<float> myDiv;
    std::vector<float> myP;

    std::vector<float> myVort;
    std::vector<float> myVortFx;
    std::vector<float> myVortFy;
};

#endif
