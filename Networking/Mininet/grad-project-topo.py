"""
Addapted from the 2 host/2 switch example provided by mininet
"""
from mininet.topo import Topo
import string

class MyTopo( Topo ):
    "Topography for Project 5"

    def __init__( self ):
        # Initialize topology
        Topo.__init__( self )

        # Add hosts and switches
        A_to_L = string.ascii_uppercase[:12]
        hosts = {}
        switches = {}

        for i, l in enumerate(A_to_L):
            hosts[l] = self.addHost( 'h' + str(i) )
            switches[l] = self.addSwitch( 's' + str(i) )
            # Every host is connected to the corresponding switch
            self.addLink( hosts[l], switches[l] )

        # Add all of the other links defined by the network. There are 19 in
        # all
        links = [
            ('A', 'B'), ('A', 'C'), ('A', 'D'), ('B', 'C'),
            ('B', 'D'), ('B', 'E'), ('C', 'D'), ('E', 'F'),
            ('E', 'G'), ('F', 'G'), ('F', 'K'), ('F', 'L'),
            ('G', 'H'), ('G', 'I'), ('G', 'J'), ('H', 'I'),
            ('H', 'J'), ('I', 'J'), ('K', 'L'),
        ]
        for l, r  in links:
            self.addLink( switches[l], switches[r] )

topos = { 'mytopo': ( lambda: MyTopo() ) }
