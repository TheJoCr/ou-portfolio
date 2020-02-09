DIR=`dirname "$0"`

# First, start the controller, using a spanning tree algorithm to prevent
# circulation of packets
cat << EOF > pox_test.mn
sh ~/pox/pox.py \
  --verbose \
  openflow.discovery \
  forwarding.l2_multi \
  openflow.spanning_tree --no-flood --hold-down \
  > pox.log 2>&1 & 
sh echo "Started POX Controller.. Waiting for completion"
sh sleep 5
sh echo "Launching Ping All"
pingall
sh echo "And validating iperf"
iperf h1 h2
EOF

# Note, during testing I found a bug in an optimization of the spanning tree
# algorithm. That is fixed locally, so this might not run elsewhere.

# Run the actual mininet network emulation, using the pox controller.  Here, we
# validate that pingall works.
sudo mn \
  --custom $DIR/grad-project-topo.py \
  --topo mytopo \
  --link tc,bw=10,delay=10ms \
  --controller=remote,ip=127.0.0.1:6633 \
  --test none \
  --post pox_test.mn

# Kill Pox.
echo "Restarting POX before validating static path insertion"
sudo pkill -9 -f '/home/mininet/pox/pox.py'
