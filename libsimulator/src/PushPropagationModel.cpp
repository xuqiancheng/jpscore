// SPDX-License-Identifier: LGPL-3.0-or-later
#include "PushPropagationModel.hpp"
#include "PushPropagationModelData.hpp"
#include "PushPropagationModelUpdate.hpp"
#include "GenericAgent.hpp"
#include "GeometricFunctions.hpp"
#include "Macros.hpp"
#include "OperationalModel.hpp"
#include "SimulationError.hpp"
#include <algorithm>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

PushPropagationModel::PushPropagationModel(double pushoutStrength, uint64_t rng_seed)
    : pushoutStrength(pushoutStrength), gen(rng_seed)
{
}

OperationalModelType PushPropagationModel::Type() const
{
    return OperationalModelType::PUSH_PROPAGATION_MODEL;
}

OperationalModelUpdate PushPropagationModel::ComputeNewPosition(
    double dT,
    const GenericAgent& ped,
    const CollisionGeometry& geometry,
    const NeighborhoodSearchType& neighborhoodSearch) const
{
    auto neighborhood = neighborhoodSearch.GetNeighboringAgents(ped.pos, _cutOffRadius);
    const auto& boundary = geometry.LineSegmentsInApproxDistanceTo(ped.pos);

    // Remove any agent from the neighborhood that is obstructed by geometry and the current
    // agent
    neighborhood.erase(
        std::remove_if(
            std::begin(neighborhood),
            std::end(neighborhood),
            [&ped, &boundary](const auto& neighbor) {
                if(ped.id == neighbor.id) {
                    return true;
                }
                const auto agent_to_neighbor = LineSegment(ped.pos, neighbor.pos);
                if(std::find_if(
                       boundary.cbegin(),
                       boundary.cend(),
                       [&agent_to_neighbor](const auto& boundary_segment) {
                           return intersects(agent_to_neighbor, boundary_segment);
                       }) != boundary.end()) {
                    return true;
                }

                return false;
            }),
        std::end(neighborhood));
 
    // calculate the push force from neighbors
    const auto pushForceHuman = std::accumulate(
        std::begin(neighborhood),
        std::end(neighborhood),
        Point{},
        [&ped, this](const auto& res, const auto& neighbor) {
            return res + PushForceHuman(ped, neighbor);
        });

    // calculate the push force from walls
    const auto boundaryForce= std::accumulate(
        boundary.cbegin(),
        boundary.cend(),
        Point(0, 0),
        [this, &ped](const auto& acc, const auto& element) {
            return acc + PushForceWall(ped, element);
        });
    printf("boundary force is (%f,%f).\n", boundaryForce.x, boundaryForce.y);

    const auto pushForce = PushForceExternal(ped) + pushForceHuman + boundaryForce;
    const auto desiredSpeed = DesiredSpeed(ped, dT);

    const auto& model = std::get<PushPropagationModelData>(ped.model);
    const auto v0 = model.velocity;
    // Get the weight of pedestrians from experiments
    double weight = 80;
    const auto velocity = desiredSpeed + (pushForce / weight) * dT;
    if (ped.id == 1)
    {
        //printf("The velocity of agent 1 at time %f is %f.\n",model.simTime, model.velocity.Norm());
    }
     auto direction = velocity.Normalized();
    if(direction == Point{}) {
        direction = ped.orientation;
    }

    // update the position, speed, and direction
    return PushPropagationModelUpdate{
        .position = ped.pos + velocity * dT, .velocity = velocity, .orientation = direction,.simTime=model.simTime+dT};
};

void PushPropagationModel::ApplyUpdate(const OperationalModelUpdate& upd, GenericAgent& agent)
    const
{
    const auto& update = std::get<PushPropagationModelUpdate>(upd);
    auto& model = std::get<PushPropagationModelData>(agent.model);
    agent.pos = update.position;
    agent.orientation = update.orientation;
    model.velocity = update.velocity;
    model.simTime = update.simTime;
}

Point PushPropagationModel::UpdateDirection(
    const GenericAgent& ped,
    const Point& calculatedDirection,
    double dt) const
{
    const auto& model = std::get<PushPropagationModelData>(ped.model);
    const Point desiredDirection = (ped.destination - ped.pos).Normalized();
    const Point actualDirection = ped.orientation;
    Point updatedDirection;

    if(desiredDirection.ScalarProduct(calculatedDirection) *
           desiredDirection.ScalarProduct(actualDirection) <
       0) {
        updatedDirection = calculatedDirection;
    } else {
        // Compute the rate of change of direction (Eq. 7)
        const Point directionDerivative =
            (calculatedDirection.Normalized() - actualDirection) / model.reactionTime;
        updatedDirection = actualDirection + directionDerivative * dt;
    }

    return updatedDirection.Normalized();
}

void PushPropagationModel::CheckModelConstraint(
    const GenericAgent& agent,
    const NeighborhoodSearchType& neighborhoodSearch,
    const CollisionGeometry& geometry) const
{
    const auto& model = std::get<PushPropagationModelData>(agent.model);
    const auto r = model.radius;
    constexpr double rMin = 0.;
    constexpr double rMax = 2.;
    validateConstraint(r, rMin, rMax, "radius", true);

    const auto strengthNeighborRepulsion = model.strengthNeighborRepulsion;
    constexpr double snMin = 0.;
    constexpr double snMax = 20.;
    validateConstraint(strengthNeighborRepulsion, snMin, snMax, "strengthNeighborRepulsion", false);

    const auto rangeNeighborRepulsion = model.rangeNeighborRepulsion;
    constexpr double rnMin = 0.;
    constexpr double rnMax = 5.;
    validateConstraint(rangeNeighborRepulsion, rnMin, rnMax, "rangeNeighborRepulsion", true);

    const auto buff = model.wallBufferDistance;
    constexpr double buffMin = 0.;
    constexpr double buffMax = 1.;
    validateConstraint(buff, buffMin, buffMax, "wallBufferDistance", false);

    const auto v0 = model.v0;
    constexpr double v0Min = 0.;
    constexpr double v0Max = 10.;
    validateConstraint(v0, v0Min, v0Max, "v0");

    const auto timeGap = model.timeGap;
    constexpr double timeGapMin = 0.;
    constexpr double timeGapMax = 10.;
    validateConstraint(timeGap, timeGapMin, timeGapMax, "timeGap", true);

    const auto anticipationTime = model.anticipationTime;
    constexpr double anticipationTimeMin = 0.0;
    constexpr double anticipationTimeMax = 5.0;
    validateConstraint(
        anticipationTime, anticipationTimeMin, anticipationTimeMax, "anticipationTime");

    const auto reactionTime = model.reactionTime;
    constexpr double reactionTimeMin = 0.0;
    constexpr double reactionTimeMax = 1.0;
    validateConstraint(reactionTime, reactionTimeMin, reactionTimeMax, "reactionTime", true);

    const auto neighbors = neighborhoodSearch.GetNeighboringAgents(agent.pos, 2);
    for(const auto& neighbor : neighbors) {
        if(agent.id == neighbor.id) {
            continue;
        }
        const auto& neighbor_model = std::get<PushPropagationModelData>(neighbor.model);
        const auto contanctdDist = r + neighbor_model.radius;
        const auto distance = (agent.pos - neighbor.pos).Norm();
        if(contanctdDist >= distance) {
            throw SimulationError(
                "Model constraint violation: Agent {} too close to agent {}: distance {}",
                agent.pos,
                neighbor.pos,
                distance);
        }
    }

    const auto lineSegments = geometry.LineSegmentsInDistanceTo(r, agent.pos);
    if(std::begin(lineSegments) != std::end(lineSegments)) {
        throw SimulationError(
            "Model constraint violation: Agent {} too close to geometry boundaries, distance "
            "<= {}",
            agent.pos,
            r);
    }
}

std::unique_ptr<OperationalModel> PushPropagationModel::Clone() const
{
    return std::make_unique<PushPropagationModel>(*this);
}

double
PushPropagationModel::OptimalSpeed(
    const GenericAgent& ped,
    double spacing,
    double time_gap) const
{
    const auto& model = std::get<PushPropagationModelData>(ped.model);
    const double min_spacing = 0.0;
    return std::min(std::max(spacing / time_gap, min_spacing), model.v0);
}

double PushPropagationModel::GetSpacing(
    const GenericAgent& ped1,
    const GenericAgent& ped2,
    const Point& direction) const
{
    const auto& model1 = std::get<PushPropagationModelData>(ped1.model);
    const auto& model2 = std::get<PushPropagationModelData>(ped2.model);
    const auto distp12 = ped2.pos - ped1.pos;
    const auto inFront = direction.ScalarProduct(distp12) >= 0;
    if(!inFront) {
        return std::numeric_limits<double>::max();
    }

    const auto left = direction.Rotate90Deg();
    const auto l = model1.radius + model2.radius;
    const bool inCorridor = std::abs(left.ScalarProduct(distp12)) <= l;
    if(!inCorridor) {
        return std::numeric_limits<double>::max();
    }
    return distp12.Norm() - l;
}

Point PushPropagationModel::CalculateInfluenceDirection(
    const Point& desiredDirection,
    const Point& predictedDirection) const
{
    // Eq. (5)
    const Point orthogonalDirection = Point(-desiredDirection.y, desiredDirection.x).Normalized();
    const double alignment = orthogonalDirection.ScalarProduct(predictedDirection);
    Point influenceDirection = orthogonalDirection;
    if(fabs(alignment) < J_EPS) {
        // Choose a random direction (left or right)
        if(gen() % 2 == 0) {
            influenceDirection = -orthogonalDirection;
        }
    } else if(alignment > 0) {
        influenceDirection = -orthogonalDirection;
    }
    return influenceDirection;
}

Point PushPropagationModel::NeighborRepulsion(
    const GenericAgent& ped1,
    const GenericAgent& ped2) const
{
    const auto& model1 = std::get<PushPropagationModelData>(ped1.model);
    const auto& model2 = std::get<PushPropagationModelData>(ped2.model);

    const auto distp12 = ped2.pos - ped1.pos;
    const auto [distance, ep12] = distp12.NormAndNormalized();
    const double adjustedDist = distance - (model1.radius + model2.radius);

    // Pedestrian movement and desired directions
    const auto& e1 = ped1.orientation;
    const auto& d1 = (ped1.destination - ped1.pos).Normalized();
    const auto& e2 = ped2.orientation;

    // Check perception range (Eq. 1)
    const auto inPerceptionRange = d1.ScalarProduct(ep12) >= 0 || e1.ScalarProduct(ep12) >= 0;
    if(!inPerceptionRange)
        return Point(0, 0);

    const double S_Gap =
        (model1.velocity - model2.velocity).ScalarProduct(ep12) * model1.anticipationTime;
    double R_dist = adjustedDist - S_Gap;
    R_dist = std::max(R_dist, 0.0); // Clamp to zero if negative

    // Interaction strength (Eq. 3 & 4)
    constexpr double alignmentBase = 1.0;
    constexpr double alignmentWeight = 0.5;
    const double alignmentFactor = alignmentBase + alignmentWeight * (1.0 - d1.ScalarProduct(e2));
    const double interactionStrength = model1.strengthNeighborRepulsion * alignmentFactor *
                                       std::exp(-R_dist / model1.rangeNeighborRepulsion);
    const auto newep12 = distp12 + model2.velocity * model2.anticipationTime; // e_ij(t+ta)

    // Compute adjusted influence direction
    const auto influenceDirection = CalculateInfluenceDirection(d1, newep12);
    return influenceDirection * interactionStrength;
}

Point PushPropagationModel::HandleWallAvoidance(
    const Point& direction,
    const Point& agentPosition,
    double agentRadius,
    const std::vector<LineSegment>& boundary,
    double wallBufferDistance) const
{
    const double criticalWallDistance = wallBufferDistance + agentRadius;
    const double influenceStartDistance =
        2.0 * criticalWallDistance; // Smoothing earlier. The constant is chosen randomly.

    auto nearestWallIt = std::min_element(
        boundary.cbegin(), boundary.cend(), [&agentPosition](const auto& wall1, const auto& wall2) {
            const auto distanceVector1 = agentPosition - wall1.ShortestPoint(agentPosition);
            const auto distanceVector2 = agentPosition - wall2.ShortestPoint(agentPosition);
            return distanceVector1.Norm() < distanceVector2.Norm();
        });

    if(nearestWallIt != boundary.end()) {
        const auto closestPoint = nearestWallIt->ShortestPoint(agentPosition);
        const auto distanceVector = agentPosition - closestPoint;
        const auto [perpendicularDistance, directionAwayFromBoundary] =
            distanceVector.NormAndNormalized();

        // Always check if too close to wall, regardless of movement direction
        if(perpendicularDistance < criticalWallDistance) {
            const auto wallVector = nearestWallIt->p2 - nearestWallIt->p1;
            const auto wallDirection = wallVector.Normalized();

            // Get parallel component of current direction
            const auto parallelComponent = wallDirection * direction.ScalarProduct(wallDirection);

            const auto newDirection =
                parallelComponent + directionAwayFromBoundary * pushoutStrength;
            return newDirection.Normalized();
        }
        // Check if within influence range
        else if(perpendicularDistance < influenceStartDistance) {
            const auto dotProduct = direction.ScalarProduct(directionAwayFromBoundary);

            // Only modify direction if moving towards wall
            if(dotProduct < 0) {
                const auto wallVector = nearestWallIt->p2 - nearestWallIt->p1;
                const auto wallDirection = wallVector.Normalized();

                if(perpendicularDistance <= criticalWallDistance) {
                    // At or closer than critical distance: enforce parallel movement
                    const auto parallelComponent =
                        wallDirection * direction.ScalarProduct(wallDirection);
                    return parallelComponent.Normalized();
                } else {
                    // Between influence start and critical distance: smooth transition
                    // Calculate influence factor: 0 at influenceStartDistance, 1 at
                    // criticalWallDistance.
                    const double influenceFactor =
                        (influenceStartDistance - perpendicularDistance) /
                        (influenceStartDistance - criticalWallDistance);

                    const auto parallelComponent =
                        wallDirection * direction.ScalarProduct(wallDirection);
                    const auto perpendicularComponent = direction - parallelComponent;

                    // Gradually reduce the perpendicular component based on distance
                    const auto newDirection =
                        parallelComponent + perpendicularComponent * (1.0 - influenceFactor);
                    return newDirection.Normalized();
                }
            }
        }
    }
    return direction;
}

Point PushPropagationModel::PushForceExternal(const GenericAgent& ped) const
{
    double force = 0;
    const auto& model = std::get<PushPropagationModelData>(ped.model);
    // todo: force should be a time series from input
    double alpha = 5;
    if (ped.id == 1 && model.simTime<0.2)
    {
        force = 250;
        /*printf(
            "external force of agent %d at time %f is %f.\n",
            ped.id,
            model.simTime,
            force);*/
    } 
    return Point(0, force * alpha);
}

Point PushPropagationModel::PushForceHuman(const GenericAgent& ped1, const GenericAgent& ped2) const
{
    const auto& model1 = std::get<PushPropagationModelData>(ped1.model);
    const auto& model2 = std::get<PushPropagationModelData>(ped2.model);

    const auto distp12 = ped2.pos - ped1.pos;
    const auto [distance, ep12] = distp12.NormAndNormalized();
    const double adjustedDist = distance - (model1.radius + model2.radius);
    // no contact no force
    if(adjustedDist >= 0)
        return Point();
   // calculate the strength of the pushing force
    double kPush = 20;
    double DPush = 0.05;
    double forceStrength = kPush * std::exp(-adjustedDist / DPush);
    return -ep12 * forceStrength;
}

Point PushPropagationModel::PushForceWall(
    const GenericAgent& ped,
    const LineSegment& boundary_segment) const
{
    const auto pt = boundary_segment.ShortestPoint(ped.pos);
    const auto dist_vec = pt - ped.pos;
    const auto [dist, e_iw] = dist_vec.NormAndNormalized();
    const auto& model = std::get<PushPropagationModelData>(ped.model);
    const auto l = model.radius;
    double adjustedDist = dist - l;
    // no contact no force
    if(adjustedDist >= 0)
        return Point();
    double kPush = 20;
    double DPush = 0.05;
    double forceStrength = kPush * std::exp(-adjustedDist / DPush);
    return -e_iw * forceStrength;
}

Point PushPropagationModel::DesiredSpeed(const GenericAgent& ped,  double dT) const
{
    const auto& model = std::get<PushPropagationModelData>(ped.model);
    const Point v0 = model.velocity;
    const Point deceDirection = -v0.Normalized();
    // constant deceleration
    double deceStrength= 2;
    Point Deceleration = deceDirection * deceStrength;
    Point desiredSpeed = v0 + Deceleration * dT;
    if (desiredSpeed.ScalarProduct(v0) < 0)
    {
        return Point();
    }
    return desiredSpeed;
}
