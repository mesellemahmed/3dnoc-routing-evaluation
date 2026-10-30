// Route-computation and selection units of the evaluated routing classes for
// an 8x8x4 mesh (3-bit X, 3-bit Y, 2-bit Z coordinates), synthesised with
// Yosys onto the NanGate 45 nm open cell library to compare their logic cost.
// Port encoding (one-hot, 7 bits): {down, up, south, north, west, east, local}.
// Each module is the per-input-port logic a router replicates for every input.

// ---------------------------------------------------------------- helpers
module dor_xyz (input [2:0] cx, cy, input [1:0] cz, input [2:0] dx, dy, input [1:0] dz,
                output reg [6:0] port);
  always @* begin
    if (dx > cx) port = 7'b0000010;
    else if (dx < cx) port = 7'b0000100;
    else if (dy > cy) port = 7'b0001000;
    else if (dy < cy) port = 7'b0010000;
    else if (dz > cz) port = 7'b0100000;
    else if (dz < cz) port = 7'b1000000;
    else port = 7'b0000001;
  end
endmodule

module dor_zxy (input [2:0] cx, cy, input [1:0] cz, input [2:0] dx, dy, input [1:0] dz,
                output reg [6:0] port);
  always @* begin
    if (dz > cz) port = 7'b0100000;
    else if (dz < cz) port = 7'b1000000;
    else if (dx > cx) port = 7'b0000010;
    else if (dx < cx) port = 7'b0000100;
    else if (dy > cy) port = 7'b0001000;
    else if (dy < cy) port = 7'b0010000;
    else port = 7'b0000001;
  end
endmodule

module dor_zyx (input [2:0] cx, cy, input [1:0] cz, input [2:0] dx, dy, input [1:0] dz,
                output reg [6:0] port);
  always @* begin
    if (dz > cz) port = 7'b0100000;
    else if (dz < cz) port = 7'b1000000;
    else if (dy > cy) port = 7'b0001000;
    else if (dy < cy) port = 7'b0010000;
    else if (dx > cx) port = 7'b0000010;
    else if (dx < cx) port = 7'b0000100;
    else port = 7'b0000001;
  end
endmodule

module minimal_mask (input [2:0] cx, cy, input [1:0] cz, input [2:0] dx, dy, input [1:0] dz,
                     output [6:0] mask);
  assign mask = {dz < cz, dz > cz, dy < cy, dy > cy, dx < cx, dx > cx,
                 (dx == cx) && (dy == cy) && (dz == cz)};
endmodule

// Index (0..6) of the largest of three 5-bit credit counts among a candidate
// mask; used as the generic credit-based selection function.
module max3_sel #(parameter W = 5) (input [6:0] mask, input [7*W-1:0] metric, input min_mode,
                                   output reg [6:0] pick);
  integer i;
  reg [W-1:0] best;
  reg found;
  always @* begin
    pick = 7'b0;
    best = min_mode ? {W{1'b1}} : {W{1'b0}};
    found = 1'b0;
    for (i = 0; i < 7; i = i + 1)
      if (mask[i] && (!found || (min_mode ? metric[i*W +: W] < best : metric[i*W +: W] > best))) begin
        best = metric[i*W +: W];
        pick = 7'b1 << i;
        found = 1'b1;
      end
  end
endmodule

// ---------------------------------------------------------------- classes
module rc_xyz (input [2:0] cx, cy, input [1:0] cz, input [2:0] dx, dy, input [1:0] dz,
               output [6:0] port);
  dor_xyz u (cx, cy, cz, dx, dy, dz, port);
endmodule

module rc_zxy (input [2:0] cx, cy, input [1:0] cz, input [2:0] dx, dy, input [1:0] dz,
               output [6:0] port);
  dor_zxy u (cx, cy, cz, dx, dy, dz, port);
endmodule

module rc_nf (input [2:0] cx, cy, input [1:0] cz, input [2:0] dx, dy, input [1:0] dz,
              input [34:0] credits, output [6:0] port);
  wire [6:0] neg = {dz < cz, 1'b0, dy < cy, 1'b0, dx < cx, 1'b0, 1'b0};
  wire [6:0] pos = {1'b0, dz > cz, 1'b0, dy > cy, 1'b0, dx > cx, 1'b0};
  wire here = (dx == cx) && (dy == cy) && (dz == cz);
  wire [6:0] adm = here ? 7'b0000001 : (|neg ? neg : pos);
  max3_sel #(5) s (adm, credits, 1'b0, port);
endmodule

module rc_dyad (input [2:0] cx, cy, input [1:0] cz, input [2:0] dx, dy, input [1:0] dz,
                input [34:0] credits, output [6:0] port);
  wire [6:0] neg = {dz < cz, 1'b0, dy < cy, 1'b0, dx < cx, 1'b0, 1'b0};
  wire [6:0] pos = {1'b0, dz > cz, 1'b0, dy > cy, 1'b0, dx > cx, 1'b0};
  wire here = (dx == cx) && (dy == cy) && (dz == cz);
  wire [6:0] adm = here ? 7'b0000001 : (|neg ? neg : pos);
  // Deterministic port: lowest set bit of the admissible set.
  wire [6:0] det = adm & (~adm + 7'd1);
  reg [4:0] det_credits;
  integer i;
  always @* begin
    det_credits = 5'd0;
    for (i = 0; i < 7; i = i + 1) if (det[i]) det_credits = credits[i*5 +: 5];
  end
  wire congested = det_credits < 5'd4;
  wire [6:0] adaptive;
  max3_sel #(5) s (adm, credits, 1'b0, adaptive);
  assign port = congested ? adaptive : det;
endmodule

// Hamiltonian (snake) labels and label-monotone admissible set; minimal moves
// preferred, credit-based choice among them.
module ham_label (input [2:0] x, y, input [1:0] z, output [7:0] lab);
  wire [5:0] idx = y[0] ? {y, ~x} : {y, x};
  wire [5:0] idx2 = z[0] ? ~idx : idx;
  assign lab = {z, idx2};
endmodule

module rc_ham (input [2:0] cx, cy, input [1:0] cz, input [2:0] dx, dy, input [1:0] dz,
               input [34:0] credits, output [6:0] port);
  wire [7:0] lc, ld, le, lw, ln, ls, lu, lo;
  ham_label a (cx, cy, cz, lc);
  ham_label b (dx, dy, dz, ld);
  ham_label e (cx + 3'd1, cy, cz, le);
  ham_label w (cx - 3'd1, cy, cz, lw);
  ham_label n (cx, cy + 3'd1, cz, ln);
  ham_label s (cx, cy - 3'd1, cz, ls);
  ham_label u (cx, cy, cz + 2'd1, lu);
  ham_label o (cx, cy, cz - 2'd1, lo);
  wire up = ld > lc;
  function ok(input [7:0] l); ok = up ? (l > lc && l <= ld) : (l < lc && l >= ld); endfunction
  wire [6:0] exists = {cz != 2'd0, cz != 2'd3, cy != 3'd0, cy != 3'd7, cx != 3'd0, cx != 3'd7, 1'b0};
  wire [6:0] adm = {ok(lo), ok(lu), ok(ls), ok(ln), ok(lw), ok(le), 1'b0} & exists;
  wire [6:0] minm;
  minimal_mask mm (cx, cy, cz, dx, dy, dz, minm);
  wire here = (dx == cx) && (dy == cy) && (dz == cz);
  // Non-minimal fallback: the admissible neighbour whose label is closest to
  // the destination label (maximum progress along the Hamiltonian path).
  function [7:0] gap(input [7:0] l); gap = (l > ld) ? l - ld : ld - l; endfunction
  wire [55:0] gaps = {gap(lo), gap(lu), gap(ls), gap(ln), gap(lw), gap(le), 8'hff};
  wire [6:0] fallback;
  max3_sel #(8) fb (adm, gaps, 1'b1, fallback);
  wire [6:0] sel;
  max3_sel #(5) sl (adm & minm, credits, 1'b0, sel);
  assign port = here ? 7'b0000001 : ((|(adm & minm)) ? sel : fallback);
endmodule

module rc_o1t (input [2:0] cx, cy, input [1:0] cz, input [2:0] dx, dy, input [1:0] dz, input order,
               output [6:0] port);
  wire [6:0] p0, p1;
  dor_xyz a (cx, cy, cz, dx, dy, dz, p0);
  dor_zyx b (cx, cy, cz, dx, dy, dz, p1);
  assign port = order ? p1 : p0;
endmodule

// Two-phase oblivious routing (ROMM / Valiant): the intermediate node travels
// in the header; the phase bit is set on arrival at the intermediate node.
module rc_twophase (input [2:0] cx, cy, input [1:0] cz, input [2:0] dx, dy, input [1:0] dz,
                    input [2:0] ix, iy, input [1:0] iz, input phase_in,
                    output [6:0] port, output phase_out);
  assign phase_out = phase_in || ((ix == cx) && (iy == cy) && (iz == cz));
  wire [2:0] tx = phase_out ? dx : ix;
  wire [2:0] ty = phase_out ? dy : iy;
  wire [1:0] tz = phase_out ? dz : iz;
  dor_xyz u (cx, cy, cz, tx, ty, tz, port);
endmodule

// Minimal fully adaptive routing with a Duato XYZ escape VC; max-credit
// selection on the adaptive VC, escape when no adaptive VC is free.
module rc_fa (input [2:0] cx, cy, input [1:0] cz, input [2:0] dx, dy, input [1:0] dz, input escaped,
              input [34:0] credits, input [6:0] adaptive_free, output [6:0] port, output use_escape);
  wire [6:0] minm, esc;
  minimal_mask mm (cx, cy, cz, dx, dy, dz, minm);
  dor_xyz d (cx, cy, cz, dx, dy, dz, esc);
  wire [6:0] cand = minm & adaptive_free;
  wire [6:0] sel;
  max3_sel #(5) s (cand, credits, 1'b0, sel);
  assign use_escape = escaped || ~|cand;
  assign port = use_escape ? esc : sel;
endmodule

// Regional congestion awareness (RCA-1D): per output direction, the
// aggregate congestion is the mean of the local metric and the value received
// from the neighbour in that direction (8-bit side-band per link).

// Q-routing selection and update datapath (16-bit Q values read from the
// per-destination table; the table itself is reported as storage).

// Thermal-aware selection: minimum 8-bit temperature reading of the
// downstream routers (side-band), escape as in fa.
module rc_ta (input [2:0] cx, cy, input [1:0] cz, input [2:0] dx, dy, input [1:0] dz, input escaped,
              input [55:0] temp, input [6:0] adaptive_free, output [6:0] port, output use_escape);
  wire [6:0] minm, esc;
  minimal_mask mm (cx, cy, cz, dx, dy, dz, minm);
  dor_xyz d (cx, cy, cz, dx, dy, dz, esc);
  wire [6:0] cand = minm & adaptive_free;
  wire [6:0] sel;
  max3_sel #(8) s (cand, temp, 1'b1, sel);
  assign use_escape = escaped || ~|cand;
  assign port = use_escape ? esc : sel;
endmodule

// Elevator-First route computation: XY towards the elevator carried in the
// header (4-bit elevator index -> coordinates), vertical at the column, XY in
// the destination layer; the virtual network follows the vertical direction.
module rc_ef (input [2:0] cx, cy, input [1:0] cz, input [2:0] dx, dy, input [1:0] dz, input [3:0] elev,
              output [6:0] port, output [1:0] vn);
  wire [2:0] ex = {elev[1:0], 1'b1};
  wire [2:0] ey = {elev[3:2], 1'b1};
  wire same = dz == cz;
  wire [2:0] tx = same ? dx : ex;
  wire [2:0] ty = same ? dy : ey;
  wire [6:0] p;
  dor_xyz u (cx, cy, cz, tx, ty, same ? cz : dz, p);
  assign port = p;
  assign vn = same ? 2'b11 : (dz > cz ? 2'b01 : 2'b10);
endmodule

// First-Last: positive-first towards the elevator, positive-last in the
// destination layer, three virtual networks; credit-based choice among the
// adaptive pair.
module rc_fl (input [2:0] cx, cy, input [1:0] cz, input [2:0] dx, dy, input [1:0] dz, input [3:0] elev,
              input [1:0] vn_in, input [34:0] credits, output [6:0] port, output [1:0] vn_out);
  wire other = dz != cz;
  wire [2:0] tx = other ? {elev[1:0], 1'b1} : dx;
  wire [2:0] ty = other ? {elev[3:2], 1'b1} : dy;
  wire ne = tx > cx, nn = ty > cy, nw = tx < cx, ns = ty < cy;
  wire at_col = !(ne || nn || nw || ns);
  wire [6:0] vert = {dz < cz, dz > cz, 5'b0};
  wire [6:0] posm = {3'b0, nn, 1'b0, ne, 1'b0};
  wire [6:0] negm = {2'b0, ns, 1'b0, nw, 2'b0};
  wire pos_first = other && (vn_in == 2'd0) && (ne || nn);
  wire neg_phase = other ? !pos_first : ((nw || ns) && vn_in != 2'd2);
  wire [6:0] adm = (dx == cx && dy == cy && dz == cz) ? 7'b0000001 :
                   (other && at_col) ? vert : (pos_first ? posm : (neg_phase ? negm : posm));
  max3_sel #(5) s (adm, credits, 1'b0, port);
  assign vn_out = (other && at_col) ? 2'd1 : (pos_first ? 2'd0 : (neg_phase ? 2'd1 : 2'd2));
endmodule

// AdEle elevator selection at the source (k = 4 candidates): EWMA costs
// (a = 0.2 ~ 13/64), relative-cost skip test rewritten without division,
// round-robin pointer and low-traffic override.
module adele_sel (input clk, input rst, input [15:0] u_rand, input fb_valid, input [1:0] fb_idx,
                  input [15:0] fb_te, input [7:0] dist0, dist1, dist2, dist3, input new_pkt,
                  output reg [1:0] choice);
  reg [15:0] c [0:3];
  reg [1:0] ptr;
  integer i;
  wire [17:0] total = c[0] + c[1] + c[2] + c[3];
  wire low = (c[0] < 16'd26) && (c[1] < 16'd26) && (c[2] < 16'd26) && (c[3] < 16'd26);
  function [1:0] argmin4(input [7:0] a, b, cc, d);
    argmin4 = (a <= b && a <= cc && a <= d) ? 2'd0 : (b <= cc && b <= d) ? 2'd1 : (cc <= d) ? 2'd2 : 2'd3;
  endfunction
  // skip_i: u * total < (4 c_i - total) * (1 - xi), with xi ~ 1/16.
  function skip(input [15:0] ci);
    reg signed [19:0] excess;
    reg [35:0] lhs, rhs;
    begin
      excess = {2'b0, ci, 2'b0} - {2'b0, total};
      lhs = u_rand * total;
      rhs = (excess > 0) ? ((excess * 15) >> 4) << 16 : 36'd0;
      skip = (excess > 0) && (lhs < rhs);
    end
  endfunction
  wire [1:0] i0 = ptr, i1 = ptr + 2'd1, i2 = ptr + 2'd2, i3 = ptr + 2'd3;
  wire [1:0] rr = !skip(c[i0]) ? i0 : !skip(c[i1]) ? i1 : !skip(c[i2]) ? i2 : i3;
  always @(posedge clk) begin
    if (rst) begin
      for (i = 0; i < 4; i = i + 1) c[i] <= 16'd0;
      ptr <= 2'd0;
      choice <= 2'd0;
    end else begin
      if (fb_valid) c[fb_idx] <= c[fb_idx] - ((c[fb_idx] * 13) >> 6) + ((fb_te * 13) >> 6);
      if (new_pkt) begin
        choice <= low ? argmin4(dist0, dist1, dist2, dist3) : rr;
        if (!low) ptr <= rr + 2'd1;
      end
    end
  end
endmodule

// ---------------------------------------------------------------- context
// One virtual-channel input buffer (8 flits of 64 bits), flip-flop based.
module vc_fifo #(parameter W = 64, D = 8) (input clk, input rst, input push, input pop,
                                           input [W-1:0] din, output [W-1:0] dout, output empty,
                                           output full);
  reg [W-1:0] mem [0:D-1];
  reg [2:0] head, tail;
  reg [3:0] count;
  assign dout = mem[head];
  assign empty = count == 4'd0;
  assign full = count == D;
  always @(posedge clk) begin
    if (rst) begin head <= 3'd0; tail <= 3'd0; count <= 4'd0; end
    else begin
      if (push && !full) begin mem[tail] <= din; tail <= tail + 3'd1; end
      if (pop && !empty) head <= head + 3'd1;
      count <= count + (push && !full) - (pop && !empty);
    end
  end
endmodule

// 7x7 crossbar for 64-bit flits (one-hot select per output).
module xbar7 (input [447:0] in_flat, input [48:0] sel_flat, output [447:0] out_flat);
  genvar o, i;
  for (o = 0; o < 7; o = o + 1) begin : outp
    wire [63:0] acc [0:7];
    assign acc[0] = 64'd0;
    for (i = 0; i < 7; i = i + 1) begin : inp
      assign acc[i+1] = acc[i] | ({64{sel_flat[o*7+i]}} & in_flat[i*64 +: 64]);
    end
    assign out_flat[o*64 +: 64] = acc[7];
  end
endmodule

// ------------------------------------------------- split shared/per-port units
// RCA-1D aggregation, one per router: six 7-bit aggregates from the local
// occupancy of each output and the value received from that neighbour.
module rca_agg (input [41:0] local_occ, input [41:0] remote_agg, output [41:0] agg_out);
  genvar g;
  for (g = 0; g < 6; g = g + 1) begin : a
    wire [7:0] s = {1'b0, local_occ[g*7 +: 7]} + {1'b0, remote_agg[g*7 +: 7]};
    assign agg_out[g*7 +: 7] = s[7:1];
  end
endmodule

module rca_select (input [2:0] cx, cy, input [1:0] cz, input [2:0] dx, dy, input [1:0] dz, input escaped,
                   input [41:0] agg, input [6:0] adaptive_free, output [6:0] port, output use_escape);
  wire [6:0] minm, esc;
  minimal_mask mm (cx, cy, cz, dx, dy, dz, minm);
  dor_xyz d (cx, cy, cz, dx, dy, dz, esc);
  wire [6:0] cand = minm & adaptive_free;
  wire [6:0] sel;
  max3_sel #(7) s ({cand}, {agg, 7'h7f}, 1'b1, sel);
  assign use_escape = escaped || ~|cand;
  assign port = use_escape ? esc : sel;
endmodule

// Q-routing: per-port selection over the 16-bit Q row of the destination and
// per-output update datapath; the Q table is counted as storage.
module qr_select (input [2:0] cx, cy, input [1:0] cz, input [2:0] dx, dy, input [1:0] dz, input escaped,
                  input [111:0] qrow, input [6:0] adaptive_free, output [6:0] port, output use_escape);
  wire [6:0] minm, esc;
  minimal_mask mm (cx, cy, cz, dx, dy, dz, minm);
  dor_xyz d (cx, cy, cz, dx, dy, dz, esc);
  wire [6:0] cand = minm & adaptive_free;
  wire [6:0] sel;
  max3_sel #(16) s (cand, qrow, 1'b1, sel);
  assign use_escape = escaped || ~|cand;
  assign port = use_escape ? esc : sel;
endmodule

module qr_update (input [15:0] q_old, input [15:0] wait_cycles, input [15:0] q_next_min,
                  output [15:0] q_new);
  wire [16:0] target = wait_cycles + q_next_min + 16'd1;
  wire [17:0] sum = q_old + target;
  assign q_new = sum[16:1];
endmodule
