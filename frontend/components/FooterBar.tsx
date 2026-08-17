import { Info, Trash2, ArrowRight, FlaskConical } from "lucide-react";

export default function FooterBar({
  onClear,
  onConfirm,
  organismSupportsMechanisms = true,
  canOptIntoMechanisms = false,
  mechanismsOptedIn = false,
  onToggleMechanismOptIn,
}: {
  onClear: () => void;
  onConfirm: () => void;
  /** Resolved from the organism's tier plus any explicit opt-in. */
  organismSupportsMechanisms?: boolean;
  /** True for Tier 4-6, where the mechanism flow is off by default. */
  canOptIntoMechanisms?: boolean;
  mechanismsOptedIn?: boolean;
  onToggleMechanismOptIn?: (next: boolean) => void;
}) {
  return (
    <div className="sticky bottom-0 z-10 border-t border-[#E5E7EB] bg-white/95 backdrop-blur px-6 py-3">
      <div className="flex flex-col gap-2 sm:flex-row sm:items-center sm:justify-between">
        <p className="flex items-center gap-2 text-[11.5px] text-slate-500">
          <Info className="h-3.5 w-3.5 shrink-0 text-slate-400" />
          {organismSupportsMechanisms
            ? "Please verify the above information before proceeding."
            : "Gene information only — mechanism analysis is not enabled for this organism."}
        </p>
        <div className="flex shrink-0 items-center gap-2">
          {/*
            Tier 4-6 are opt-in rather than blocked. The rulebooks describe
            mammalian RNA biology, so enabling them for a plant, virus or
            bacterium is the user's call to make explicitly — but the option
            is right here rather than hidden.
          */}
          {canOptIntoMechanisms && (
            <label className="flex cursor-pointer items-center gap-1.5 rounded border border-amber-300 bg-amber-50 px-2.5 py-1.5 text-[11.5px] font-medium text-amber-800 transition-colors hover:bg-amber-100">
              <input
                type="checkbox"
                checked={mechanismsOptedIn}
                onChange={(e) => onToggleMechanismOptIn?.(e.target.checked)}
                className="h-3 w-3 accent-amber-600"
              />
              <FlaskConical className="h-3 w-3" />
              Enable mechanism analysis anyway
            </label>
          )}
          <button
            onClick={onClear}
            className="flex items-center gap-1.5 rounded border border-slate-300 px-3 py-1.5 text-[12px] font-medium text-slate-600 hover:bg-slate-50 transition-all duration-200 hover:border-slate-400 hover:shadow-sm active:translate-y-0"
          >
            <Trash2 className="h-3 w-3" />
            Clear All
          </button>
          {organismSupportsMechanisms && (
            <button
              onClick={onConfirm}
              className="flex items-center gap-1.5 rounded bg-brand px-4 py-1.5 text-[12px] font-medium text-white hover:bg-brand-dark transition-all duration-200 hover:shadow-md hover:shadow-brand/20 hover:-translate-y-0.5 active:translate-y-0"
            >
              Confirm &amp; Proceed
              <ArrowRight className="h-3 w-3" />
            </button>
          )}
        </div>
      </div>
    </div>
  );
}
