// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once

#include "Point.hpp"

struct PushPropagationModelData {
    double strengthNeighborRepulsion{};
    double rangeNeighborRepulsion{};
    double wallBufferDistance{0.1}; // buff distance of agent to wall
    double anticipationTime{0.5}; // anticipation time
    double reactionTime{0.1}; // reaction time to update direction
    Point velocity{};
    double timeGap{1.06};
    double v0{1.2};
    double radius{0.15};
    double simTime{0}; // total simulation time
};

template <>
struct fmt::formatter<PushPropagationModelData> {

    constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

    template <typename FormatContext>
    auto format(const PushPropagationModelData& m, FormatContext& ctx) const
    {
        return fmt::format_to(
            ctx.out(),
            "PushPropagationModel[strengthNeighborRepulsion={}, "
            "rangeNeighborRepulsion={}, wallBufferDistance={}, "
            "timeGap={}, v0={}, radius={}, reactionTime={}, anticipationTime={}, velocity={}, simTime={}])",
            m.strengthNeighborRepulsion,
            m.rangeNeighborRepulsion,
            m.wallBufferDistance,
            m.timeGap,
            m.v0,
            m.radius,
            m.reactionTime,
            m.anticipationTime,
            m.velocity,
            m.simTime);
    }
};
