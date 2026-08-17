"use client";

import { Info } from "lucide-react";
import {
  getOrganism,
  needsOrganismBanner,
  organismBannerMessage,
} from "@/lib/organisms";

/**
 * Advisory banner for any organism outside Tier 1.
 *
 * Tier 1 gets nothing — there is no caveat to make for a clinical species.
 * Tier 2-3 run the full flow but the rulebooks and delivery-precedent tables
 * are built on human evidence. Tier 4-6 need an explicit opt-in because the
 * rulebooks describe mammalian RNA biology.
 *
 * Nothing is hidden on the strength of this banner; it sets expectations.
 */
export default function OrganismCapabilityBanner({
  organismId,
}: {
  organismId: string;
}) {
  if (!needsOrganismBanner(organismId)) return null;

  const organism = getOrganism(organismId);
  const message = organismBannerMessage(organismId);
  if (!organism || !message) return null;

  const optIn = organism.tier >= 4;
  const tone = optIn
    ? "border-amber-200 bg-amber-50 text-amber-900"
    : "border-sky-200 bg-sky-50 text-sky-900";

  return (
    <div
      role="note"
      className={`flex items-start gap-2 rounded border px-3 py-2 text-[11.5px] leading-snug ${tone}`}
    >
      <Info className="mt-0.5 h-3.5 w-3.5 shrink-0 opacity-70" />
      <p>{message}</p>
    </div>
  );
}
