===================
Taking Calibrations
===================

This document describes the PypeIt-recommended way
to take spectroscopic calibrations for instruments
supported by this DRP.

LRIS
----

For each detector and binning take a set of 11+ bias frames.


Blue
++++

.. _lris_blue_arc:

Arcs
::::

Use *all* lamps and integrate for 1s after waiting for 60s for
the lamps to warm-up.

Flats
:::::

Best:

  - Illumination: Obtain a series of ~3 exposures through the mask on the bright
    twilight sky.  Aim for at least 1000 counts above bias.
  - Pixel:  Obtain a series of ~5 exposures with the mask *removed*
    on the twilight sky.  Aim for ~20,000 counts per exposure while
    keeping integration times <20s each.  Dither by 1' between
    exposures

Possibly acceptable:

  - Take a series (~5) of dome flats for as long as you can afford
    without saturating the red portions of the spectrum (<40,000).

Not acceptable:

  - Internal flats


Red
+++

Arcs
::::

Same as :ref:`lris_blue_arc`.

Flats
:::::

Best:

  - Take a series (~5) dome flats aiming for a
    peak of ~20,000 counts.

Avoid internal flats at all costs.