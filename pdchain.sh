_pdchain_usage() { cat << EOF
usage: pdchain [CMD]
Shortcut for running tools in the PDchain workspace.
Replace [CMD] with any PDchain command line tools.

For example:
pdchain rfd_chain 100-100 --outprefix outputs/example0 --relax
pdchain mpros --help

CMD special options:
    activate                Activates default PDchain environment within the
                            current shell. Allows direct input commands without
                            pdchain prefix.

                            E.g. pdchain activate
                            rfd_chain 100-100 --outprefix outputs/example1

                            Keep in mind that the PDchain environment has GNU
                            coreutils, so commands such as ls, cd, mkdir, etc.
                            will default to them. This is true for bash, sed,
                            grep, awk, and pymol as well.

    shell                   Enters subshell with PDchain environment. Allows
                            direct input commands just like with "activate".
                            Because you're in a subshell, non-exported functions
                            and variables will not be inherited. You can exit
                            the subshell with the "exit" command.

                            WARNING: Do not put this in your ~/.bashrc! From my
                            experience (as of March 2026), it gets trapped in
                            some infinite loop. Might be a bug with pixi.

    list                    List available command line tools from PDchain.

    --help                  Show this message

If you want to run RFdiffusion, LigandMPNN, etc. directly, you must specify the
environment name with -e, then specify the command with the same name.
E.g. pdchain -e rfdiffusion rfdiffusion [...]
E.g. pdchain -e ligandmpnn ligandmpnn [...]
EOF
}

_pdchain_list() {
    echo "List of available command line tools:"
    echo "esmfold_relax evo_rfd_chain mpnn_score mpros omegafold_relax rfd_chain smiles_to_params af3webtools bashscript calc_cst_rmsd calc_fixedres_rmsd calc_scafscore clean_pdb dss_disicl filter_pool gen_atompair_csts gen_contigs gen_coord_csts gen_invrots gen_natbiasjson pdb_to_fasta track_lineage update_contigs hqa_convert" | tr ' ' '\n' | sort | xargs -n4 printf "    %-19s %-19s %-19s %-s\n"
}

pdchain() {
    _pdchain_toml="/Users/johncheng/Workspaces/PDchain/pixi.toml"

    if [ "$1" = "activate" ]; then
        eval "$(pixi shell-hook -m "$_pdchain_toml")"
    elif [ "$1" = "shell" ]; then
        pixi shell -m "$_pdchain_toml"
    elif [ "$1" = "list" ]; then
        _pdchain_list
    elif [ "$1" = "--help" ] || [ "$#" -eq 0 ]; then
        _pdchain_usage
    else
        pixi run -m "$_pdchain_toml" "$@"
    fi
}

if [ -n "$BASH_VERSION" ] ; then
    export -f pdchain _pdchain_usage _pdchain_list
fi
