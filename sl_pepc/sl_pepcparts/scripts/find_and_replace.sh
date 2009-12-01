#!/bin/sh

#  
#  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
#  
#  file: scripts/find_and_replace.sh
#  timestamp: 2009-10-22 09:31:12 +0200
#  

#set -xv

# find-and-replace script ...

# ... with sed
sed_delim="&"  # sed delimiter for find-and-replace (should not occur in file/directory paths)
sed_ww_begin="\\\\<"  # control sequence for matching only whole words
sed_ww_end="\\\\>"

far_script_sed_begin()
{
  printf ""
}

far_script_sed_end()
{
  printf ""
}

far_script_sed_mod_list()
{
  local ww=$1
  local list=$2
  local mod_in=$3
  local mod_out=$4

  local _ifs=$IFS
  IFS=":"
  set $mod_in
  local in_prefix=$1
  local in_suffix=$2
  set $mod_out
  local out_prefix=$1
  local out_suffix=$2
  IFS=$_ifs

  if [ "${ww}" != "0" ] ; then
    local ww_begin=${sed_ww_begin}
    local ww_end=${sed_ww_end}
  fi

  exec < ${list}

  while read name ; do
    [ -n "$name" ] && printf "s${sed_delim}${ww_begin}${in_prefix}${name}${in_suffix}${ww_end}${sed_delim}${out_prefix}${name}${out_suffix}${sed_delim}g\n"
  done
}

far_script_sed_single()
{
  local ww=$1
  local find=$2
  local replace=$3

  if [ "${ww}" != "0" ] ; then
    local ww_begin=${sed_ww_begin}
    local ww_end=${sed_ww_end}
  fi

  printf "s${sed_delim}${ww_begin}${find}${ww_end}${sed_delim}${replace}${sed_delim}g\n"
}


# ... with perl
perl_delim="!"
perl_ww_begin="\\\\b"
perl_ww_end="\\\\b"

far_script_perl_begin()
{
  printf "while (<>) {\n"
}

far_script_perl_end()
{
  printf "print; }\n"
}

far_script_perl_mod_list()
{
  local ww=$1
  local list=$2
  local mod_in=$3
  local mod_out=$4

  local _ifs=$IFS
  IFS=":"
  set $mod_in
  local in_prefix=$1
  local in_suffix=$2
  set $mod_out
  local out_prefix=$1
  local out_suffix=$2
  IFS=$_ifs

  if [ "${ww}" != "0" ] ; then
    local ww_begin=${perl_ww_begin}
    local ww_end=${perl_ww_end}
  fi

  exec < ${list}

  while read name ; do
    [ -n "$name" ] && printf "s${perl_delim}${ww_begin}${in_prefix}${name}${in_suffix}${ww_end}${perl_delim}${out_prefix}${name}${out_suffix}${perl_delim}g;\n"
  done
}

far_script_perl_single()
{
  local ww=$1
  local find=$2
  local replace=$3

  if [ "${ww}" != "0" ] ; then
    local ww_begin=${perl_ww_begin}
    local ww_end=${perl_ww_end}
  fi

  printf "s${perl_delim}${ww_begin}${find}${ww_end}${perl_delim}${replace}${perl_delim}g;\n"
}

# FIXME: sed \< \> ability testen!

in_cfg=$1
in_cmd=$2
shift 2

[ -z "${in_cfg}" ] && in_cfg=sed

case ${in_cfg} in
  sed)
    cmd_exec="sed -f"
    cmd_begin=far_script_sed_begin
    cmd_end=far_script_sed_end
    cmd_mod_list=far_script_sed_mod_list
    cmd_single=far_script_sed_single
    ;;
  perl)
    cmd_exec="perl"
    cmd_begin=far_script_perl_begin
    cmd_end=far_script_perl_end
    cmd_mod_list=far_script_perl_mod_list
    cmd_single=far_script_perl_single
    ;;
  *)
    printf "error unknown option: '${in_cfg}'\n" >&2
    exit
    ;;
esac


case ${in_cmd} in
  exec)
    ${cmd_exec} $*
    ;;
  begin)
    ${cmd_begin} $*
    ;;
  end)
    ${cmd_end} $*
    ;;
  mod_list)
    ${cmd_mod_list} $*
    ;;
  single)
    ${cmd_single} $*
    ;;
  *)
    printf "error unknown command: '${in_cmd}'\n" >&2
    exit
    ;;
esac
