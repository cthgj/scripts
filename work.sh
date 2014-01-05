#!/usr/bin/env bash

xterm -e 'top; bash' &
xterm -e 'free -h; bash' &  && wmctrl -r 'xterm' -b add,maximized_vert &
# xterm -e '/bin/bash --rcfile ~/3.rc' &
# xterm -e '/bin/bash --rcfile ~/4.rc' &
