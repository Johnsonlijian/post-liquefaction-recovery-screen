# R39 Module 3 — validation + decision payload

## (a) Leave-one-coordinate-out generalisation of the recovery rate
full within-coordinate rate: 1.95 qc1N/month
LOO rates: min 1.15, max 2.87, all positive = True
-> recovery direction is not driven by any single coordinate.

## (b) Tolerable next-event PGA as the ground recovers (median loose site, M6 aftershock)
  critical layer z=3.0 m, water 1.5 m, start qc1Ncs=37.5, beta=1.95/mo
   t=   0 mo  qc1Ncs=  37.5  tolerable PGA (FS=1) = 0.110 g
   t=   3 mo  qc1Ncs=  43.4  tolerable PGA (FS=1) = 0.116 g
   t=   6 mo  qc1Ncs=  49.2  tolerable PGA (FS=1) = 0.122 g
   t=  12 mo  qc1Ncs=  60.9  tolerable PGA (FS=1) = 0.135 g
   t=  24 mo  qc1Ncs=  84.3  tolerable PGA (FS=1) = 0.170 g
   t=  60 mo  qc1Ncs= 154.5  tolerable PGA (FS=1) = 0.592 g
   (linear beta cannot continue indefinitely: it would imply implausibly dense states;
    physically recovery decelerates as log-time aging -> another reason the form matters.)
-> One year of recovery raises tolerable PGA from 0.110 g to 0.135 g (+26 mg): recovery buys little demand resistance per year.
   Decision use: the clock tells you the demand level below which rebuilt reserve suffices,
   and the (long) time needed to tolerate a given aftershock PGA.

## (c) Manifestation cross-check (Feb 2011), honest limitations
  pairs with mapped manifestation polygon: 6/35
  fraction mapped as liquefaction (primary): 0.17
  NOTE: 'no_polygon' is not verified absence; mapping is incomplete, so this is
  a supportive consistency check, not a calibrated out-of-sample test. The clock
  predicts deficit-driven re-liquefaction at the Feb-2011 demand for ~97% of sites,
  consistent with the observed widespread manifestation.