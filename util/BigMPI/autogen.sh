#! /bin/sh
autoreconf=${AUTORECONF:-autoreconf}
$autoreconf ${autoreconf_args:-"-vif"}
