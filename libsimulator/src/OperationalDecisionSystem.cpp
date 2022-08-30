#include "OperationalDecisionSystem.hpp"
#include "OperationalModel.hpp"
#include <memory>

OperationalDecisionSystem::OperationalDecisionSystem(std::unique_ptr<OperationalModel>&& model)
    : _model(std::move(model))
{
}

void OperationalDecisionSystem::Run(
    double dT,
    double /*t_in_sec*/,
    const NeighborhoodSearch& neighborhoodSearch,
    const CollisionGeometry& geometry,
    std::vector<std::unique_ptr<Pedestrian>>& agents) const
{

    std::vector<std::optional<PedestrianUpdate>> updates(agents.size(), std::nullopt);
    std::transform(
        agents.begin(),
        agents.end(),
        updates.begin(),
        [this, dT, &neighborhoodSearch, &geometry](
            auto&& agent) -> std::optional<PedestrianUpdate> {
            return _model->ComputeNewPosition(dT, *agent, geometry, neighborhoodSearch);
        });

    auto agent_iter = agents.begin();
    std::for_each(updates.begin(), updates.end(), [this, &agent_iter](auto&& update) {
        if(update) {
            _model->ApplyUpdate(*update, **agent_iter);
        }
        ++agent_iter;
    });
}
