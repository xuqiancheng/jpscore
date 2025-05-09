# SPDX-License-Identifier: LGPL-3.0-or-later

from dataclasses import dataclass
from random import randint

import jupedsim.native as py_jps


@dataclass(kw_only=True)
class PushPropagationModel:
    """Push Propagation Model (PPM).

    Attributes:
        pushout_strength: The pushout mechanism ensures agents maintain a safe distance
        from walls by adding a small outward component to their movement when within the
        critical wall distance. This outward component, scaled by `pushoutStrength`,
        combines with the parallel component of the agent's direction to create smooth,
        gliding behavior along walls.
        rng_seed: seed value of internally used rng. If not explicitly set this
            value will be chosen randomly.
    """

    pushout_strength: float = 0.3
    rng_seed: int = randint(0, 2**64 - 1)


@dataclass(kw_only=True)
class PushPropagationModelAgentParameters:
    """
    Agent parameters for Push Propagation Model (PPM).

    .. note::

        Instances of this type are copied when creating the agent, you can safely
        create one instance of this type and modify it between calls to `add_agent`

        E.g.:

        .. code:: python

            positions = [...] # List of initial agent positions
            params = AnticipationVelocityModelAgentParameters(desired_speed=0.9) # all agents are slower
            for p in positions:
                params.position = p
                sim.add_agent(params)

    Attributes:
        position: Position of the agent.
        time_gap: Time constant that describe how fast pedestrian close gaps.
        desired_speed: Maximum speed of the agent.
        radius: Radius of the agent.
        journey_id: Id of the journey the agent follows.
        stage_id: Id of the stage the agent targets.
        strength_neighbor_repulsion: Strength of the repulsion from neighbors
        range_neighbor_repulsion: Range of the repulsion from neighbors
        wall_buffer_distance: Buffer distance of agents to the walls.
        anticipation_time: Anticipation time of an agent.
        reaction_time: reaction time of an agent to change its direction.
        external_force: external force from the punching bag.
    """

    position: tuple[float, float] = (0.0, 0.0)
    time_gap: float = 1.06
    desired_speed: float = 1.2
    radius: float = 0.2
    journey_id: int = 0
    stage_id: int = 0
    strength_neighbor_repulsion: float = 8.0
    range_neighbor_repulsion: float = 0.1
    wall_buffer_distance: float = 0.1
    anticipation_time: float = 1.0
    reaction_time: float = 0.3
    external_force: float=0
    weight: float=80
    height: float=1.8
    forceAmplify: float=1

    def as_native(
        self,
    ) -> py_jps.PushPropagationModelAgentParameters:
        return py_jps.PushPropagationModelAgentParameters(
            position=self.position,
            time_gap=self.time_gap,
            desired_speed=self.desired_speed,
            radius=self.radius,
            journey_id=self.journey_id,
            stage_id=self.stage_id,
            strength_neighbor_repulsion=self.strength_neighbor_repulsion,
            range_neighbor_repulsion=self.range_neighbor_repulsion,
            wall_buffer_distance=self.wall_buffer_distance,
            anticipation_time=self.anticipation_time,
            reaction_time=self.reaction_time,
            external_force=self.external_force,
            weight=self.weight,
            height=self.height,
            forceAmplify=self.forceAmplify,
        )


class PushPropagationModelState:
    def __init__(self, backing):
        self._obj = backing

    @property
    def time_gap(self) -> float:
        return self._obj.time_gap

    @time_gap.setter
    def time_gap(self, time_gap):
        self._obj.time_gap = time_gap

    @property
    def desired_speed(self) -> float:
        """desired Speed of this agent."""
        return self._obj.desired_speed

    @desired_speed.setter
    def desired_speed(self, desired_speed):
        self._obj.desired_speed = desired_speed

    @property
    def radius(self) -> float:
        """Radius of this agent."""
        return self._obj.radius

    @radius.setter
    def radius(self, radius):
        self._obj.radius = radius

    @property
    def strength_neighbor_repulsion(self) -> float:
        """Strength of the repulsion from neighbors of this agent."""
        return self._obj.strength_neighbor_repulsion

    @strength_neighbor_repulsion.setter
    def strength_neighbor_repulsion(self, strength_neighbor_repulsion):
        self._obj.strength_neighbor_repulsion = strength_neighbor_repulsion

    @property
    def range_neighbor_repulsion(self) -> float:
        """Range of the repulsion from neighbors of this agent."""
        return self._obj.range_neighbor_repulsion

    @range_neighbor_repulsion.setter
    def range_neighbor_repulsion(self, range_neighbor_repulsion):
        self._obj.range_neighbor_repulsion = range_neighbor_repulsion

    @property
    def wall_buffer_distance(self):
        """Wall buffer distance of agent to walls."""
        return self._obj.wall_buffer_distance

    @wall_buffer_distance.setter
    def wall_buffer_distance(self, wall_buffer_distance):
        self._obj.wall_buffer_distance = wall_buffer_distance

    @property
    def anticipation_time(self) -> float:
        """Anticipation time of this agent."""
        return self._obj.anticipation_time

    @anticipation_time.setter
    def anticipation_time(self, anticipation_time):
        self._obj.anticipation_time = anticipation_time

    @property
    def reaction_time(self) -> float:
        """Reaction time of this agent."""
        return self._obj.reaction_time

    @reaction_time.setter
    def reaction_time(self, reaction_time):
        self._obj.reaction_time = reaction_time

    @property
    def external_force(self) -> float:
        """external force of this agent."""
        return self._obj.external_force

    @external_force.setter
    def external_force(self, external_force):
        self._obj.external_force = external_force

        
    @property
    def weight(self) -> float:
        """weight of this agent."""
        return self._obj.weight

    @weight.setter
    def weight(self, weight):
        self._obj.weight = weight

    @property
    def height(self) -> float:
        """height of this agent."""
        return self._obj.height

    @height.setter
    def height(self, height):
        self._obj.height = height

    @property
    def forceAmplify(self) -> float:
        """forceAmplify of this agent."""
        return self._obj.forceAmplify

    @forceAmplify.setter
    def forceAmplify(self, forceAmplify):
        self._obj.forceAmplify = forceAmplify