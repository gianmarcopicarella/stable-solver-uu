//
// Created by Gianmarco Picarella on 11/04/23.
//

#ifndef GAMEPHYSICSPROJECT_VISUALIZER_H
#define GAMEPHYSICSPROJECT_VISUALIZER_H

class StableSolver;

class Visualizer {
public:
    explicit Visualizer(const StableSolver  * aSolver) :
            mySolver(aSolver) {}

    virtual void operator()() = 0;

protected:
    const StableSolver * mySolver{};
};


class VelocityVisualizer : public Visualizer
{
public:
    explicit VelocityVisualizer(const StableSolver  * aSolver) : Visualizer(aSolver) {}
    void operator()() override;
};

class DensityVisualizer : public Visualizer
{
public:
    explicit DensityVisualizer(const StableSolver  * aSolver) : Visualizer(aSolver) {}
    void operator()() override;
};


#endif //GAMEPHYSICSPROJECT_VISUALIZER_H
