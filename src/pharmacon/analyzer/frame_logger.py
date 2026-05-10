"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

Module for logging information about frames in a Molecular Dynamics (MD) trajectory.

This module provides a context manager utility to log details of trajectory frames during
analysis in `MDAnalysis.Universe` objects. Users can control the range and frequency
of logging using slicing parameters and a configurable interval.

This is useful for debugging, monitoring, or tracking specific trajectory processing steps.
"""


import logging
from MDAnalysis import Universe
from contextlib import contextmanager




__all__ = [
    "log_every_frame",
]


@contextmanager
def log_every_frame(u: Universe,
                    logger: logging.Logger,
                    *,
                    start: int | None = None,
                    stop: int | None = None,
                    step: int = 1,
                    every: int = 1,
                    prefix: str = "frame") -> None:
    """
    Context manager for logging each accessed frame of a trajectory in the MDAnalysis
    universe. This utility is useful for efficiently tracking and reporting the progress or
    status of trajectory processing at regular intervals. Optionally restricts the logging
    to a specific slice of frames within the trajectory and supports specifying logging
    granularity.

    :param u: Universe to track, encapsulating the trajectory data being analyzed.
    :type u: Universe
    :param logger: Logger object to emit debug-level messages about processed frames.
    :type logger: logging.Logger
    :param start: Specifies the first frame in the trajectory to consider. If `None`, starts from
                  the beginning of the trajectory. Defaults to `None`.
    :type start: int | None
    :param stop: Specifies the frame index immediately after the last frame to consider. If `None`,
                 processes until the end of the trajectory. Defaults to `None`.
    :type stop: int | None
    :param step: Interval defining how many frames to skip between each processed frame.
                 Must be positive and defaults to 1.
    :type step: int
    :param every: Frequency for logging frames that meet the slice criteria during processing. For example,
                  with `every=10`, a log message will be emitted only for every 10th processed frame.
                  Defaults to 1.
    :type every: int
    :param prefix: Prefix string to prepend to log messages for identification of logged frames.
                   Defaults to "frame".
    :type prefix: str
    :return: None
    """

    # Preserve existing transformations
    prev = list(getattr(u.trajectory, "_transformations", []))

    # Normalize slice parameters
    start0 = 0 if start is None else int(start)
    stop0 = None if stop is None else int(stop)
    step0 = 1 if step in (None, 0) else int(step)
    every0 = max(1, int(every))

    def in_slice(f: int) -> bool:
        """
        Determines whether a given integer falls within a specific slice defined
        by the variables `start0`, `stop0`, and `step0`.

        The function checks if the input integer `f` lies within the range
        [start0, stop0) (inclusive start, exclusive stop) and validates
        whether it adheres to the step interval constraints specified by
        `step0`, if applicable.

        :param f: The integer to check against the slice constraints.
        :type f: int
        :return: A boolean value indicating whether the integer `f` satisfies
                 all slice conditions.
        :rtype: bool
        """
        if f < start0:
            return False
        if stop0 is not None and f >= stop0:
            return False
        if step0 and ((f - start0) % step0 != 0):
            return False
        return True

    state = {
        "prev_f": None,          # last *touched* frame
        "prev_t": None,          # last touched time
        "count_in_slice": 0,     # number of completed frames in slice
        "last_done_f": None,     # last frame that was ACTUALLY processed
        "last_done_t": None,     # time of that frame
    }

    def _tap(ts):
        """
        Processes a timestamp object and updates the internal state for frame tracking,
        logging information as required based on configured parameters.

        :param ts: The timestamp object to process. It is expected to have attributes:
                   `frame` (an integer representing the current frame) and `time`
                   (a floating-point time associated with the frame, which may be None).
        :return: Returns the input `ts` after processing.
        """
        f = int(getattr(ts, "frame", 0))
        t = getattr(ts, "time", None)
        t = None if t is None else float(t)

        # First frame ever seen
        if state["prev_f"] is None:
            state["prev_f"] = f
            state["prev_t"] = t
            return ts

        # A new frame was accessed → previous frame is DONE
        done_f = state["prev_f"]
        done_t = state["prev_t"]

        if in_slice(done_f):
            n = state["count_in_slice"]
            if n % every0 == 0:
                if done_t is None:
                    logger.debug("%s #%d", prefix, done_f)
                else:
                    logger.debug("%s #%d (t=%.3f ps)", prefix, done_f, done_t)

            state["count_in_slice"] += 1
            state["last_done_f"] = done_f
            state["last_done_t"] = done_t

        # Update previous → current
        state["prev_f"] = f
        state["prev_t"] = t
        return ts

    try:
        u.trajectory._transformations = prev + [_tap]
        yield
    finally:
        # Log ONLY the last ACTUALLY processed frame
        f = state["last_done_f"]
        t = state["last_done_t"]

        if f is not None:
            if t is None:
                logger.debug("%s #%d DONE", prefix, f)
            else:
                logger.debug("%s #%d DONE (t=%.3f ps)", prefix, f, t)

        # Restore original transformations
        u.trajectory._transformations = prev
