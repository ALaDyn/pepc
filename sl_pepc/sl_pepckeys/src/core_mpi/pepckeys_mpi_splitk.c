/*
 *  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
 *  
 *  file: src/core_mpi/pepckeys_mpi_splitk.c
 *  timestamp: 2009-11-13 18:17:04 +0100
 *  
 */

/* Planloses isend/irecv kann zu deadlock führen WENN mpi-implementation bei wait[all|any] nicht alle requests prüft sondern nur EINEN
   -> HAHA bei wait müssen (!?!) alle Verbindungen überprüft werden

    p1, p2, p3 -> p0
    p4, p5, p6 -> p7
    p0 -> p7
    p7 -> p0

    - alle senden jeweils genau x Elemente und haben auch nur DIE x Elemente frei
    - p0 und p7 empfangen jeweils zuerst von [p1|p2|p3] und [p4|p5|p6] und sind voll
    -> deadlock zwischen p0 und p7, sind voll und müssen aber vom jeweiligen gegenüber empfangen

  Favorit: alltoall
   * jeder Prozess empfängt so viel wie er versendet, deadlockfrei möglich wenn genügend Platz zum empfangen vorhanden
   * Platz zum empfangen schaffen, indem auf send_buffer verteilt wird
   * komplette Verteilung der send_buffer -> ein Prozess muss mehr empfangen als er versendet
   * komplette Verteilung der send_buffer reduzieren?
   * Reduzierung darf nicht direkt/indirekt in deadlock enden

   * 1. beliebige Reduzierung bis Platz ausreicht
   Gegenbeispiel:
     - p0 <- p1, p2 <- p3, p0 <-> p2
     - alle Prozesse habe x Elemente frei und wollen x Elemente senden
     - Reduzierung (beliebig) auf p0 <- p1, p2 <- p3
     - im nächsten Schritt wieder: p0 <- p1, p2 <- p3, p0 <-> p2
     - p0 und p2 sind voll, deadlock

   * 2. Reduzierung um Zyklen
     - Beweis: bei 1 Element send_buffer pro Zielprozess, ist immer mindestens ein Zyklus enthalten
     - zuerst doppelte (hin und zurück) Kanten zwischen Knoten auflösen
     - 0 Kanten entfernen -> Netzwerk enstanden (gerichteter Graph ohne Mehrfachkanten)
     - so lange es noch starke Zusammenhangskomponenten gibt, in einer starken Zusammenhangskomponenten Zyklen suchen und eliminieren
     -> 1 Knoten verdoppeln (Quelle und Senke bilden), maximalen Fluss bestimmen (Ford und Fulkerson) und entfernen

 Implementierung:
  - MPI_Allgather besser als MPI_Gather + MPI_Scatter
  
  - Zeiten von k2c_equal vs. k2c_part = 0.84 : 1 (ohne Daten) bzw. 0.95 : 1(mit Daten)
  -> "aufwendige" key2class_part Funktion benötigt 12% von ewalk (bei 2^22 elements/node)
  -> "einfache" pepckeys_key2class_equal Funktion benötigt 8% von ewalk (bei 2^22 elements/node)
  
 Erweiterung:
  - "fertige" Prozesse aus Kommunikator (und damit auch aus Berechnung) entfernen
  - sinnvoll wenn durch ungleiche Verteilung Prozesse unterschiedlich schnell fertig werden
  - leicht möglich nach sample_select?

 Rückgabewert:
  x == 0 -> success
  x < 0 -> beliebiger Fehler (code x)
  x > 0 -> deadlock, die erst x Elemente sind korrekt bearbeitet worden

 Rückgabeparameter:
  nne = neue Anzahl Elemente in s und sa
  nue = Anzahl nicht bearbeiteter Elemente (besser nfe = Anzahl bearbeiteter Elemente?)
  
  JUMP: bei size = 32 liegt flow_reduce (28 ms) über flow_gather (14-20 ms)

*/


#include "sl_common.h"


slint pepckeys_flow_create(slint *m, slint n, slint *complete, slint *reduced, int *send_counts, int *recv_counts, int rank) /* sl_func pepckeys_flow_create */
{
  slint i, j, k, y, z;
  
  *complete = *reduced = 0;

  /* reduce multiple-edges */
  for (j = 0; j < n; ++j)
  for (i = 0; i < j; ++i)
  {
    y = m[j * n + i];
    z = m[i * n + j];
    *complete += y + z;
    
    /* compute the direct exchange between node i and j */
    k = xmin(y, z);

    /* reduce the flow by the direct exchange between node i and j */
    m[j * n + i] -= k;
    m[i * n + j] -= k;
    *reduced += 2 * k;
    
    /* update the send_/recv_counts */
    if (rank == j) send_counts[i] = recv_counts[i] = k;
    if (rank == i) send_counts[j] = recv_counts[j] = k;
  }

  send_counts[rank] = recv_counts[rank] = 0;

  return 0;
}


#define white -1
#define black -2
#define init_flow INT_MAX

slint pepckeys_dfs_visit_rbe(slint v, slint maxflow, slint *pred, slint *m, slint n, slint *reduced, int *send_counts, int *recv_counts, int rank) /* sl_func pepckeys_dfs_visit_rbe */
{
  slint i, j, k;

/*  printf("visiting %d\n", v);*/

  for (i = 0; i < n; i++)
  if ((k = m[v * n + i]) > 0)
  {
    if (pred[i] == white) /* not visited, yet */
    {
      pred[i] = v;  /* mark as grey */

      if (pepckeys_dfs_visit_rbe(i, xmin(maxflow, k), pred, m, n, reduced, send_counts, recv_counts, rank)) return 1;

      pred[i] = black;  /* mark as black */

    } else if (pred[i] != black) /* not white, not black -> it's grey, it's a back edge */
    {
      k = xmin(maxflow, k);

/*      printf("%d elements cycle:", k);*/

      j = i;
      do 
      {
        m[v * n + j] -= k;
        *reduced += k;

        if (rank == v) send_counts[j] += k;
        if (rank == j) recv_counts[v] += k;

/*        printf(" %d -> %d ", v, j);*/
        
        j = v;
        v = pred[j];

      } while (j != i);
      
/*      printf("\n");*/

      return 1;
    }
  }

/*  printf("visiting %d done\n", v);*/
  
  return 0;
}


/* komplizierter, aber schneller möglich?

 - merken, welche Kante auf dem bisherigen Weg die minimalste war, bzw. von welchem Knoten diese minimalste Kante aus genommen wurde
 - dieses FRÜHESTE ("MIT GIBS NICH", nur bei lt wird neues Minimum genommen) "minimale" Kante gilt im jeder Rekusionstiefe für den zurückliegenden aktuellen Pfad
 - wenn Rückwärtskante gefunden, minimalen Pfad entfernen ABER nicht neu starten
 - weitermachen, wenn Minimum immer noch > 0, sonst: Rekursion bis zur frühsten 0 (alten & neues frühstes Minimum) zurückgehen und ab da mit neuen Kanten weitermachen
 - unabhängig davon: wenn Startknoten i, dann nur noch Kanten zu Konten j >= i untersuchen
 - bsp: nachdem Startknoten 0 fertig ist, ist er an keinen Zyklen mehr beteiligt + muss auch nichtmehr untersucht werden
 
 -> wie wirkt sich dies auf Komplexität (Theorie) und Laufzeit (Praxis) aus?
*/

slint pepckeys_cycle_reduce_back_edges(slint *m, slint n, slint *reduced, int *send_counts, int *recv_counts, int rank) /* sl_func pepckeys_cycle_reduce_back_edges */
{
  slint i, j, pred[n], found_back_edge = 1;
  
/*  printf("reducing back edges\n");*/

  while (found_back_edge)
  {
    for (i = 0; i < n; i++) pred[i] = white;

    found_back_edge = 0;
    for (i = 0; i < n; i++)
    if (pred[i] == white)
    {
      pred[i] = i;

      if (pepckeys_dfs_visit_rbe(i, init_flow, pred, m, n, reduced, send_counts, recv_counts, rank))
      {
        found_back_edge = 1;
        break;
      }
      
      pred[i] = black;
    }
  }

  return 0;
}


slint pepckeys_flow_reduce_unbalanced(slint *m, slint n, slint *reduced, slint *max_unbalance, int *send_counts, int *recv_counts, int rank) /* sl_func pepckeys_flow_reduce_unbalanced */
{
  slint i, j, k;
  
  for (i = 0; i < n; ++i)
  for (j = 0; j < n; ++j)
  if (i != j)
  {
    k = xmin(max_unbalance[j], m[i * n + j]);
    
    if (rank == i) send_counts[j] += k;
    if (rank == j) recv_counts[i] += k;
    
    max_unbalance[j] -= k;
    m[i * n + j] -= k;
    
    *reduced += k;
  }
  
  return 0;
}


/* sx -> send_buffer, sa -> auxiliary space, nne -> final number of elements (may be less, equal or greater than s->size), nue -> final number of unfinished elements */
slint pepckeys_mpi_splitk(elements_t *s, k2c_func k2c, void *ci, elements_t *sx, elements_t *sa, slint *nne, slint *nue, int size, int rank, MPI_Comm comm) /* pepckeys_sl_proto, sl_func pepckeys_mpi_splitk */
{
  slint i, j, k, l, t;

  slint unfinished, complete, reduced, remaining;

  slint local_sb_counts[size], all_sb_counts[size * size];

  elements_t sb[size], sb_current[size];
  elements_t src, dst, end, dst_prev;

  elements_t sa_dst, sa_end;
  
  int send_counts[size], send_displs[size], recv_counts[size], recv_displs[size];
  
  slint add_free_space[size];

  MPI_Comm *current_comm = &comm;

  /* return if somethings missing */  
  if (s == NULL || ci == NULL || sx == NULL) return -1;

  /* need send_buffers with at least one element per foreign process */
  if (sx->size < size - 1) return -2;

  rti_treset(rti_tid_mpi_splitk_exact_loop_walk);
  rti_treset(rti_tid_mpi_splitk_exact_loop_flow);
  rti_treset(rti_tid_mpi_splitk_exact_loop_flow_gather);
  rti_treset(rti_tid_mpi_splitk_exact_loop_flow_create);
  rti_treset(rti_tid_mpi_splitk_exact_loop_flow_reduce);
  rti_treset(rti_tid_mpi_splitk_exact_loop_dist);
  rti_treset(rti_tid_mpi_splitk_exact_loop_dist_pre);
  rti_treset(rti_tid_mpi_splitk_exact_loop_dist_a2av);

  rti_tstart(rti_tid_mpi_splitk_exact);
  rti_tstart(rti_tid_mpi_splitk_exact_init);

  /* initials */
  j = sx->size;
  k = size - 1;
  for (i = 0; i < size; ++i)
  {
    /* init the local send_buffer counters */
    local_sb_counts[i] = 0;

    /* prepare the send_buffers */
    if (i != rank)
    {
      elem_assign_at(sx, sx->size - j, &sb[i]);
      sb[i].size = (j / k) + (j % k != 0);
      j -= sb[i].size;
      --k;

    } else elem_null(&sb[i]);
    elem_assign(&sb[i], &sb_current[i]);
  }

  /* optional: use cloned communicator, MPI_Comm_dub */

  elem_assign(s, &src);
  elem_assign(s, &dst);
  elem_assign_at(s, s->size, &end);
  
  if (sa != NULL)
  {
    elem_assign(sa, &sa_dst);
    elem_assign_at(sa, sa->size, &sa_end);

  } else
  {
    elem_null(&sa_dst);
    elem_null(&sa_end);
  }

  rti_tstop(rti_tid_mpi_splitk_exact_init);
  rti_tstart(rti_tid_mpi_splitk_exact_loop);

  /* the target-process of the current element is unkown, used to call the key2class function only ONCE per element */
  t = -1;

  /* the main loop */
  while (1)
  {
    rti_tstart(rti_tid_mpi_splitk_exact_loop_walk);

    /* distribute the elements to the send_buffer, as long as possible (elements left and target send_buffer not full) */
    while (src.keys != end.keys)
    {
      /* if the target-process of the current element is unknown, compute it */
      if (t < 0) t = (k2c)(src.keys, src.keys - s->keys, ci);

      /* is the local process the target? */
      if (t == rank)
      {
        /* if necessary, move the element on the local process */
        if (src.keys != dst.keys) elem_copy(&src, &dst);

        /* update the dst-position */
        elem_inc(&dst);
        
      } else /* the target is another process (need to send the element) */
      {
        /* break, if the according send_buffer is full */
        if (local_sb_counts[t] >= sb[t].size) break;

        /* copy the element to the according send_buffer */
        elem_copy(&src, &sb_current[t]);
        elem_inc(&sb_current[t]);
        ++local_sb_counts[t];
      }

      /* go to the next element, signal the yet unkown target-process */
      elem_inc(&src);
      t = -1;
    }

    rti_tstop(rti_tid_mpi_splitk_exact_loop_walk);
    rti_tstart(rti_tid_mpi_splitk_exact_loop_flow);

    /* cumulate the local send-buffer counts */
    local_sb_counts[rank] = 0;
    for (i = 0; i < size; ++i)
    if (i != rank) local_sb_counts[rank] += local_sb_counts[i];

    /* if we're receiving in the ordinary space */
    if (dst.keys < end.keys) local_sb_counts[rank] = (src.keys - dst.keys) - local_sb_counts[rank];  /* offer the additionally free ordinary space (beside the free space for the regular deadlock-free cycling exchange) */
    else local_sb_counts[rank] = sa_end.keys - sa_dst.keys;                                          /* offer the auxiliary space */
    
    rti_tstart(rti_tid_mpi_splitk_exact_loop_flow_gather);

    /* distribute the send_buffer counts */
    MPI_Allgather(local_sb_counts, size, sl_int_type_mpi, all_sb_counts, size, sl_int_type_mpi, comm);

    rti_tstop(rti_tid_mpi_splitk_exact_loop_flow_gather);
    
    /* extract and remove the sizes of the local auxiliary space (hidden in the local send_buffer counts) */
    for (i = 0; i < size; i++)
    {
      /* store the sizes of local auxiliary space (offered by the processes) */
      add_free_space[i] = all_sb_counts[i * size + i];
      /* remove these "hidden" informations from the send_buffer counts */
      all_sb_counts[i * size + i] = 0;
    }

    rti_tstart(rti_tid_mpi_splitk_exact_loop_flow_create);

    /* reduce the direct xchanges between two nodes */
    pepckeys_flow_create(all_sb_counts, size, &complete, &reduced, send_counts, recv_counts, rank);

    rti_tstop(rti_tid_mpi_splitk_exact_loop_flow_create);

    /* if the flow isn't empty */
    if (reduced < complete)
    {
      rti_tstart(rti_tid_mpi_splitk_exact_loop_flow_reduce);

      /* reduce cycles, simple: making dfs and reducing back edges as long as possible */
      pepckeys_cycle_reduce_back_edges(all_sb_counts, size, &reduced, send_counts, recv_counts, rank);

      rti_tstop(rti_tid_mpi_splitk_exact_loop_flow_reduce);
    }

    /* if the flow is still not empty */
    if (reduced < complete)
    {
      rti_tstart(rti_tid_mpi_splitk_exact_loop_flow_unbalance);
    
      /* detect unbalanced sends (to additionally free space) */
      pepckeys_flow_reduce_unbalanced(all_sb_counts, size, &reduced, add_free_space, send_counts, recv_counts, rank);

      rti_tstop(rti_tid_mpi_splitk_exact_loop_flow_unbalance);
    }

    rti_tstop(rti_tid_mpi_splitk_exact_loop_flow);

    /* if there are no more unfinished elements, we're done */
    if (complete == 0)
    {
      /* set the 'end' behind the last valid element in the ordinary space */
      elem_assign(&dst, &end);
      break;
    }

    /* if nothing was reduced, we have a deadlock */
    if (!reduced)
    {
      /* write back all the elements currently residing in the send_buffer */
      elem_assign(&dst, &dst_prev);
      for (i = 0; i < size; ++i)
      if (i != rank)
      {
        elem_ncopy(&sb[i], &dst_prev, local_sb_counts[i]);
        elem_add(&dst_prev, local_sb_counts[i]);
      }

      /* close a possible hole between the previously back-written elements and the original elements at 'src' */
      k = xmin(src.keys - dst_prev.keys, end.keys - src.keys);
      i = (end.keys - src.keys) - k;
      elem_add(&src, k);
      elem_ncopy(&src, &dst_prev, k);

      /* set the 'end' behind the last valid element in the ordinary space */
      elem_assign_at(&dst_prev, k, &end);

      break;
    }

    /* distribute (go with the flow!) */
    rti_tstart(rti_tid_mpi_splitk_exact_loop_dist);
    rti_tstart(rti_tid_mpi_splitk_exact_loop_dist_pre);

    /* compute displacements by sending the last x elemens of each send_buffer */
    j = k = 0;
    for (i = 0; i < size; i++)
    {
      /* send the last elements of the send_buffer */

      local_sb_counts[i] -= send_counts[i];
      elem_sub(&sb_current[i], send_counts[i]);

      send_displs[i] = j + local_sb_counts[i];
      j += sb[i].size;

      recv_displs[i] = k;
      k += recv_counts[i];
    }

    /* if there's empty ordinary space */
    if (dst.keys >= s->keys && dst.keys < end.keys)
    {
      /* prepare receiving to the ordinary space */
      elem_assign(&dst, &dst_prev);
      elem_add(&dst, k);

    } else
    {
      /* prepare receiving to the auxiliary space */
      elem_assign(&sa_dst, &dst_prev);
      elem_add(&sa_dst, k);
    }

    rti_tstop(rti_tid_mpi_splitk_exact_loop_dist_pre);
    rti_tstart(rti_tid_mpi_splitk_exact_loop_dist_a2av);

    /* distribute the keys */
#define xelem_call \
    MPI_Alltoallv(xelem_buf(sx), send_counts, send_displs, xelem_mpi_datatype, xelem_buf(&dst_prev), recv_counts, recv_displs, xelem_mpi_datatype, comm);
#include "sl_xelem_call.h"

    rti_tstop(rti_tid_mpi_splitk_exact_loop_dist_a2av);
    rti_tstop(rti_tid_mpi_splitk_exact_loop_dist);

    /* optional: remove finished processes from the cloned communicator, MPI_Comm_split, MPI_Comm_free */
    /* notice removed processes and update an original_rank to reduced_rank mapping */
  }

  rti_tstop(rti_tid_mpi_splitk_exact_loop);
  rti_tstop(rti_tid_mpi_splitk_exact);
  
  /*  */
  if (nne != NULL) *nne = (end.keys - s->keys) + (sa_dst.keys - sa->keys);
  
  if (nue != NULL) *nue = end.keys - dst.keys;

  if (complete) printf("%d here: WARNING: nothing reduced -> DEADLOCK\n", rank);

  return (complete != 0);
}
