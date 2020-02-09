DIR=`dirname "$0"`

cat << 'EOF' > setup_Controller.sh
# First, we get link names that we are going to need. Here, we go from B ( h1 -
# 00:00:00:00:00:02 ) to F ( h5 - 00:00:00:00:00:06) along switches B, E, and F
# - that is s1, s4, and s5

# We use s1-eth5 s3-eth3 and s5-eth1 to get to our destination.  In order to
# get the interface numbers (for the output) we run:
out1=`sudo ovs-ofctl show s1 | grep 'eth5' | grep -o '[0-9]' | head -n 1`
out2=`sudo ovs-ofctl show s4 | grep 'eth3' | grep -o '[0-9]' | head -n 1`
out3=`sudo ovs-ofctl show s5 | grep 'eth1' | grep -o '[0-9]' | head -n 1`

# Coming back we use s5-eth2 s4-eth2 and s1-eth1 
back1=`sudo ovs-ofctl show s5 | grep 'eth2' | grep -o '[0-9]' | head -n 1`
back2=`sudo ovs-ofctl show s4 | grep 'eth2' | grep -o '[0-9]' | head -n 1`
back3=`sudo ovs-ofctl show s1 | grep 'eth1' | grep -o '[0-9]' | head -n 1`

src=00:00:00:00:00:02
dest=00:00:00:00:00:06
# With that defined, it's easy to insert the flow table (making sure to add
# both directions)
sudo ovs-ofctl add-flow s1 dl_src=$src,dl_dst=$dest,actions=output:$out1
sudo ovs-ofctl add-flow s4 dl_src=$src,dl_dst=$dest,actions=output:$out2
sudo ovs-ofctl add-flow s5 dl_src=$src,dl_dst=$dest,actions=output:$out3

sudo ovs-ofctl add-flow s5 dl_src=$dest,dl_dst=$src,actions=output:$back1
sudo ovs-ofctl add-flow s4 dl_src=$dest,dl_dst=$src,actions=output:$back2
sudo ovs-ofctl add-flow s1 dl_src=$dest,dl_dst=$src,actions=output:$back3

# Demonstrate that it is set up. 
echo "Flow table on s1:"
sudo ovs-ofctl dump-flows s1
echo "Flow table on s4:"
sudo ovs-ofctl dump-flows s4
echo "Flow table on s5:"
sudo ovs-ofctl dump-flows s5 
EOF

cat << EOF > pox_test2.mn
sh $DIR/setup_Controller.sh
h1 ping h5
EOF

sudo mn \
  --custom $DIR/grad-project-topo.py \
  --topo mytopo \
  --controller=remote,ip=127.0.0.1:6633 \
  --test none \
  --post pox_test2.mn \
  --mac --arp

