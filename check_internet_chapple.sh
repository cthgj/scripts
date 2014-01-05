#!/bin/bash

while [ 1 ]; do
    echo "===============================" >> /home/terdon/net_status
    date >> /home/terdon/net_status
    ping -c 1 google.com 2>> /home/terdon/net_status |tail -2 | head -1 >> /home/terdon/net_status
    echo "" >> /home/terdon/net_status
    sleep 30;
done
