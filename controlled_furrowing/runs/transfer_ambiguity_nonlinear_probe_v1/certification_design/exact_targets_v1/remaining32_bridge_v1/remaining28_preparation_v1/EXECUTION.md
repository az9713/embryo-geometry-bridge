# Wave-aware execution adapter

The adapter supports the seven frozen profiles covering exactly 28 later targets. It adds no new science kernels. Existing selector and component functions keep the same trial, cache, route, binary R, exact drive and physical observer archive. Each new raw launch has a distinct wave directory. The adapter exposes complete typed raw metadata to later independent hull and global-reduction readers. It does not start those separate proof phases.

## Before actual work

The independent cold source gate is execution_preparation/independent_prelaunch/review.json. Use PASS_INDEPENDENT_REMAINING28_WAVE_EXECUTION_COLD, nonempty literal-true checks, permitted_waves wave_02 through wave_08, maximum_owned_scientific_processes=4 and actual_release_required=true. Bind source_sha256, protocol_sha256, preflight_sha256, profiles_sha256 and joined_bridge_input_admission_sha256 from the frozen preparation.

After the independent first-wave global replay genuinely finishes, root writes that replay's actual_replay/os_exit_observation.json with PASS_ROOT_OBSERVED_FIRST_WAVE_INDEPENDENT_GLOBAL_EXIT, positive literal pid and unified_exec_session, observed_exit_code=0 as a literal integer, terminal_sha256 and receipt_sha256. The new admission/later_waves_release.json uses PASS_ROOT_RELEASE_REMAINING28_AFTER_FIRST_WAVE_GLOBAL_REPLAY, the five cold binding hashes, source_review_sha256, first_wave_independent_observed_exit_sha256 and first_wave_independent_receipt_sha256. It has literal-true checks, case_order of the first wave, following_waves from the profile, selector_invocations_only_under_this_release=true and per_wave_plan_review_and_raw_admission_required=true. An unresolved first-wave certificate is retained as such. The release still requires a complete independently reviewed first-wave replay.

## Commands

Run wave_execution.py selector-commands for a declared --wave with the exact --protocol-sha256, --preflight-sha256, --review, --review-sha256, --release and --release-sha256. It prints four ordered argv lists for the existing selector and launches nothing. Run them with genuine process handles and verify the independent actual selector replays before any raw work. Schedule selectors and later replay processes within the same overall four-scientific-process budget. The raw supervisor enforces its own shared OS lock across all seven wave IDs.

For raw work, also supply --wave-admission-sha256. The canonical waves/[wave]/admission/raw_wave_admission.json status is PASS_INDEPENDENT_LATER_WAVE_SPATIAL_ADMISSION. It repeats the release bindings and SHA, includes wave_id, copied raw_worker_targets as targets, and four plans rows in canonical order. Plan rows use the existing plan/terminal/observer/review hashes plus the observed independent plan replay SHA. Set independent_component_hull_required=true and other_waves_authorized=false. Independent source review and actual release alone cannot launch this phase.

Run the raw subcommand with the same options. The native-raw entry point verifies actual process ownership and the parent's start hash. One shared lock protects launch through native child observation. The frozen scheduler tracks each armed task, drains futures through console interrupts, and reaps owned subprocesses on failure before releasing the lock. An unclosed prior supervisor blocks later raw waves. Failed or partial amounts are preserved without automatic retry.

## Later proof input

After genuinely observing the raw supervisor, root pins the exact 715-file closure in waves/[wave]/raw/root_terminal_admission.json. Use PASS_ROOT_OBSERVED_LATER_WAVE_RAW_SPATIAL_TERMINALS, wave_id, pid, unified_exec_session, observed_exit_code=0, raw_wave_admission_sha256 and sha256. The open_raw_phase API returns metadata and file references after checking all eight native child exits and matching every raw component hash reference. It does not read component bytes or arrays or grant proof status. Reuse the existing independent local hull_math.local and signed global arithmetic with these explicit case/phase arguments and their own approval. No narrowing of a hull or acceptance of a disjoint row is authorized.

The initial synthetic fixture exposed a float alias on the external release PID before any actual work. The draft is retained. The corrected fixtures reject that alias, mixed first-wave cases and sibling plan proofs. Four real inert child processes exited under the same ownership primitives. No actual selector, spatial integral, hull replay or global reduction ran.
