# NOT Observing checklist and useful scripts

## General Precautions

- Do not drink tap water at the service room
- There's a safety laser which is placed at the exit staircase of the dome. When you get out of the dome make sure to step over it (and not through it) or the safety alarm will go off
- When you’re driving up/down make sure your car’s headlights are on parking mode otherwise you’ll pollute light on the MAGIC and other telescopes (potentially destroying their detectors). Also make sure to pull down the blinders in your room in the Residencia/hostel

## Observing checklist

- Start FIES calibs (bias, halogen lamp flats) before going to dinner
- ALFOSC can be done in the morning — the script is clever enough to check which filters/grisms that had been used during the night. Note: Sky flats won’t be done this way of course, you have to do the manually.

## Good to know
- StanCam is a camera that’s used for capturing the aquisition image for placing the fiber on a target when using FIES. StanCam can also be used separately on its own.

## Making Finding Charts
There are two finding chart generators here. I use `finding_chart.py`. I added a method that can fetch photometry, deredden colors for dust extinction using [`dustmaps`](https://dustmaps.readthedocs.io/)


The `FCgenerator.py` for generating finding charts is from Johan Fynbo.


