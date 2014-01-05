#!/bin/bash

while [ 1 ]; do
    echo "===============================" >> /home/cchapple/net_status
    date >> /home/cchapple/net_status
    ping -c 1 google.com 2>> /home/cchapple/net_status |tail -2 | head -1 >> /home/cchapple/net_status
    echo "" >> /home/cchapple/net_status
    sleep 30;
done
