#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: Pi4Dqpsk Rx
# GNU Radio version: 3.10.12.0

from gnuradio import blocks
from gnuradio import digital
from gnuradio import filter
from gnuradio.filter import firdes
from gnuradio import gr
from gnuradio.fft import window
import sys
import signal
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation
from gnuradio import network
import cmath
import osmosdr
import time
import threading




class pi4dqpsk_rx(gr.top_block):

    def __init__(self):
        gr.top_block.__init__(self, "Pi4Dqpsk Rx", catch_exceptions=True)
        self.flowgraph_started = threading.Event()

        ##################################################
        # Variables
        ##################################################
        self.sps = sps = 2
        self.nfilts = nfilts = 32
        self.constel = constel = digital.constellation_dqpsk().base()
        self.constel.set_npwr(1.0)
        self.constel.gen_soft_dec_lut(8)
        self.variable_adaptive_algorithm_0 = variable_adaptive_algorithm_0 = digital.adaptive_algorithm_cma( constel, 10e-3, 1).base()
        self.samp_rate = samp_rate = 1024000
        self.rrc_taps = rrc_taps = firdes.root_raised_cosine(nfilts, nfilts, 1.0/float(sps), 0.35, 11*sps*nfilts)
        self.decim = decim = 32
        self.channel_rate = channel_rate = 36000
        self.arity = arity = 4

        ##################################################
        # Blocks
        ##################################################

        self.rtlsdr_source = osmosdr.source(
            args="numchan=" + str(1) + " " + ''
        )
        self.rtlsdr_source.set_time_unknown_pps(osmosdr.time_spec_t())
        self.rtlsdr_source.set_sample_rate(samp_rate)
        self.rtlsdr_source.set_center_freq((420.150*1e6-150*1e3), 0)
        self.rtlsdr_source.set_freq_corr(0, 0)
        self.rtlsdr_source.set_dc_offset_mode(0, 0)
        self.rtlsdr_source.set_iq_balance_mode(0, 0)
        self.rtlsdr_source.set_gain_mode(True, 0)
        self.rtlsdr_source.set_gain(15, 0)
        self.rtlsdr_source.set_if_gain(0, 0)
        self.rtlsdr_source.set_bb_gain(0, 0)
        self.rtlsdr_source.set_antenna('', 0)
        self.rtlsdr_source.set_bandwidth(0, 0)
        self.network_udp_sink_0 = network.udp_sink(gr.sizeof_char, 1, '127.0.0.1', 42000, 0, 508, False)
        self.mmse_resampler_xx_0_0 = filter.mmse_resampler_cc(0, (float(samp_rate)/(float(decim)*float(channel_rate))))
        self.freq_xlating_fir_filter_xxx_0 = filter.freq_xlating_fir_filter_ccc(decim, firdes.low_pass(1,samp_rate,12500,12500*0.2), (150*1e3), samp_rate)
        self.digital_pfb_clock_sync_xxx_0 = digital.pfb_clock_sync_ccf(sps, (2*cmath.pi/100.0), rrc_taps, nfilts, (nfilts/2), 1.5, sps)
        self.digital_map_bb_0 = digital.map_bb(constel.pre_diff_code())
        self.digital_linear_equalizer_0 = digital.linear_equalizer(15, sps, variable_adaptive_algorithm_0, True, [ ], 'corr_est')
        self.digital_fll_band_edge_cc_0 = digital.fll_band_edge_cc(sps, 0.35, 45, (cmath.pi/100.0))
        self.digital_diff_phasor_cc_0 = digital.diff_phasor_cc()
        self.digital_constellation_decoder_cb_0 = digital.constellation_decoder_cb(constel)
        self.blocks_unpack_k_bits_bb_0 = blocks.unpack_k_bits_bb(constel.bits_per_symbol())


        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_unpack_k_bits_bb_0, 0), (self.network_udp_sink_0, 0))
        self.connect((self.digital_constellation_decoder_cb_0, 0), (self.digital_map_bb_0, 0))
        self.connect((self.digital_diff_phasor_cc_0, 0), (self.digital_constellation_decoder_cb_0, 0))
        self.connect((self.digital_fll_band_edge_cc_0, 0), (self.digital_pfb_clock_sync_xxx_0, 0))
        self.connect((self.digital_linear_equalizer_0, 0), (self.digital_diff_phasor_cc_0, 0))
        self.connect((self.digital_map_bb_0, 0), (self.blocks_unpack_k_bits_bb_0, 0))
        self.connect((self.digital_pfb_clock_sync_xxx_0, 0), (self.digital_linear_equalizer_0, 0))
        self.connect((self.freq_xlating_fir_filter_xxx_0, 0), (self.mmse_resampler_xx_0_0, 0))
        self.connect((self.mmse_resampler_xx_0_0, 0), (self.digital_fll_band_edge_cc_0, 0))
        self.connect((self.rtlsdr_source, 0), (self.freq_xlating_fir_filter_xxx_0, 0))


    def get_sps(self):
        return self.sps

    def set_sps(self, sps):
        self.sps = sps
        self.set_rrc_taps(firdes.root_raised_cosine(self.nfilts, self.nfilts, 1.0/float(self.sps), 0.35, 11*self.sps*self.nfilts))

    def get_nfilts(self):
        return self.nfilts

    def set_nfilts(self, nfilts):
        self.nfilts = nfilts
        self.set_rrc_taps(firdes.root_raised_cosine(self.nfilts, self.nfilts, 1.0/float(self.sps), 0.35, 11*self.sps*self.nfilts))

    def get_constel(self):
        return self.constel

    def set_constel(self, constel):
        self.constel = constel
        self.digital_constellation_decoder_cb_0.set_constellation(self.constel)

    def get_variable_adaptive_algorithm_0(self):
        return self.variable_adaptive_algorithm_0

    def set_variable_adaptive_algorithm_0(self, variable_adaptive_algorithm_0):
        self.variable_adaptive_algorithm_0 = variable_adaptive_algorithm_0

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.freq_xlating_fir_filter_xxx_0.set_taps(firdes.low_pass(1,self.samp_rate,12500,12500*0.2))
        self.mmse_resampler_xx_0_0.set_resamp_ratio((float(self.samp_rate)/(float(self.decim)*float(self.channel_rate))))
        self.rtlsdr_source.set_sample_rate(self.samp_rate)

    def get_rrc_taps(self):
        return self.rrc_taps

    def set_rrc_taps(self, rrc_taps):
        self.rrc_taps = rrc_taps
        self.digital_pfb_clock_sync_xxx_0.update_taps(self.rrc_taps)

    def get_decim(self):
        return self.decim

    def set_decim(self, decim):
        self.decim = decim
        self.mmse_resampler_xx_0_0.set_resamp_ratio((float(self.samp_rate)/(float(self.decim)*float(self.channel_rate))))

    def get_channel_rate(self):
        return self.channel_rate

    def set_channel_rate(self, channel_rate):
        self.channel_rate = channel_rate
        self.mmse_resampler_xx_0_0.set_resamp_ratio((float(self.samp_rate)/(float(self.decim)*float(self.channel_rate))))

    def get_arity(self):
        return self.arity

    def set_arity(self, arity):
        self.arity = arity




def main(top_block_cls=pi4dqpsk_rx, options=None):
    tb = top_block_cls()

    def sig_handler(sig=None, frame=None):
        tb.stop()
        tb.wait()

        sys.exit(0)

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    tb.start()
    tb.flowgraph_started.set()

    tb.wait()


if __name__ == '__main__':
    main()
