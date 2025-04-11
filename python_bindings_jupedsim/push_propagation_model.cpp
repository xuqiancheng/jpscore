// SPDX-License-Identifier: LGPL-3.0-or-later
#include "conversion.hpp"
#include "jupedsim/push_propagation_model.h"
#include "wrapper.hpp"

#include <jupedsim/jupedsim.h>

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void init_push_propagation_model(py::module_& m)
{
    py::class_<JPS_PushPropagationModelAgentParameters>(
        m, "PushPropagationModelAgentParameters")
        .def(
            py::init([](std::tuple<double, double> position,
                        double time_gap,
                        double desired_speed,
                        double radius,
                        JPS_JourneyId journey_id,
                        JPS_StageId stage_id,
                        double strengthNeighborRepulsion,
                        double rangeNeighborRepulsion,
                        double wallBufferDistance,
                        double anticipationTime,
                        double reactionTime) {
                return JPS_PushPropagationModelAgentParameters{
                    intoJPS_Point(position),
                    journey_id,
                    stage_id,
                    time_gap,
                    desired_speed,
                    radius,
                    strengthNeighborRepulsion,
                    rangeNeighborRepulsion,
                    wallBufferDistance,
                    anticipationTime,
                    reactionTime};
            }),
            py::kw_only(),
            py::arg("position"),
            py::arg("time_gap"),
            py::arg("desired_speed"),
            py::arg("radius"),
            py::arg("journey_id"),
            py::arg("stage_id"),
            py::arg("strength_neighbor_repulsion"),
            py::arg("range_neighbor_repulsion"),
            py::arg("wall_buffer_distance"),
            py::arg("anticipation_time"),
            py::arg("reaction_time"))
        .def("__repr__", [](const JPS_PushPropagationModelAgentParameters& p) {
            return fmt::format(
                "position: {}, journey_id: {}, stage_id: {}, "
                "time_gap: {}, desired_speed: {}, radius: {}",
                "strength_neighbor_repulsion: {}, range_neighbor_repulsion: {}"
                "wall_buffer_distance: {}"
                "anticipation_time: {}, reaction_time: {}",
                intoTuple(p.position),
                p.journeyId,
                p.stageId,
                p.time_gap,
                p.v0,
                p.radius,
                p.strengthNeighborRepulsion,
                p.rangeNeighborRepulsion,
                p.wallBufferDistance,
                p.anticipationTime,
                p.reactionTime);
        });
    py::class_<JPS_PushPropagationModelBuilder_Wrapper>(m, "PushPropagationModelBuilder")
        .def(
            py::init([](double pushoutStrength, uint64_t rng_seed) {
                return std::make_unique<JPS_PushPropagationModelBuilder_Wrapper>(
                    JPS_PushPropagationModelBuilder_Create(pushoutStrength, rng_seed));
            }),
            py::kw_only(),
            py::arg("pushout_strength"),
            py::arg("rng_seed"))
        .def("build", [](JPS_PushPropagationModelBuilder_Wrapper& w) {
            JPS_ErrorMessage errorMsg{};
            auto result = JPS_PushPropagationModelBuilder_Build(w.handle, &errorMsg);
            if(result) {
                return std::make_unique<JPS_OperationalModel_Wrapper>(result);
            }
            auto msg = std::string(JPS_ErrorMessage_GetMessage(errorMsg));
            JPS_ErrorMessage_Free(errorMsg);
            throw std::runtime_error{msg};
        });
    py::class_<JPS_PushPropagationModelState_Wrapper>(m, "PushPropagationModelState")
        .def_property(
            "time_gap",
            [](const JPS_PushPropagationModelState_Wrapper& w) {
                return JPS_PushPropagationModelState_GetTimeGap(w.handle);
            },
            [](JPS_PushPropagationModelState_Wrapper& w, double time_gap) {
                JPS_PushPropagationModelState_SetTimeGap(w.handle, time_gap);
            })
        .def_property(
            "desired_speed",
            [](const JPS_PushPropagationModelState_Wrapper& w) {
                return JPS_PushPropagationModelState_GetV0(w.handle);
            },
            [](JPS_PushPropagationModelState_Wrapper& w, double desiredSpeed) {
                JPS_PushPropagationModelState_SetV0(w.handle, desiredSpeed);
            })
        .def_property(
            "radius",
            [](const JPS_PushPropagationModelState_Wrapper& w) {
                return JPS_PushPropagationModelState_GetRadius(w.handle);
            },
            [](JPS_PushPropagationModelState_Wrapper& w, double radius) {
                JPS_PushPropagationModelState_SetRadius(w.handle, radius);
            })
        .def_property(
            "strength_neighbor_repulsion",
            [](const JPS_PushPropagationModelState_Wrapper& w) {
                return JPS_PushPropagationModelState_GetStrengthNeighborRepulsion(w.handle);
            },
            [](JPS_PushPropagationModelState_Wrapper& w, double strengthNeighborRepulsion) {
                JPS_PushPropagationModelState_SetStrengthNeighborRepulsion(
                    w.handle, strengthNeighborRepulsion);
            })
        .def_property(
            "range_neighbor_repulsion",
            [](const JPS_PushPropagationModelState_Wrapper& w) {
                return JPS_PushPropagationModelState_GetRangeNeighborRepulsion(w.handle);
            },
            [](JPS_PushPropagationModelState_Wrapper& w, double rangeNeighborRepulsion) {
                JPS_PushPropagationModelState_SetRangeNeighborRepulsion(
                    w.handle, rangeNeighborRepulsion);
            })
        .def_property(
            "wall_buffer_distance",
            [](const JPS_PushPropagationModelState_Wrapper& w) {
                return JPS_PushPropagationModelState_GetWallBufferDistance(w.handle);
            },
            [](JPS_PushPropagationModelState_Wrapper& w, double wallBufferDistance) {
                JPS_PushPropagationModelState_SetWallBufferDistance(
                    w.handle, wallBufferDistance);
            })
        .def_property(
            "anticipation_time",
            [](const JPS_PushPropagationModelState_Wrapper& w) {
                return JPS_PushPropagationModelState_GetAnticipationTime(w.handle);
            },
            [](JPS_PushPropagationModelState_Wrapper& w, double anticipationTime) {
                JPS_PushPropagationModelState_SetAnticipationTime(w.handle, anticipationTime);
            })
        .def_property(
            "reaction_time",
            [](const JPS_PushPropagationModelState_Wrapper& w) {
                return JPS_PushPropagationModelState_GetReactionTime(w.handle);
            },
            [](JPS_PushPropagationModelState_Wrapper& w, double reactionTime) {
                JPS_PushPropagationModelState_SetReactionTime(w.handle, reactionTime);
            });
}
