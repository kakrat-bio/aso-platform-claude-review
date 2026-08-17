"""
SQLAlchemy models for users, profiles, saved work, and activity.
"""

from __future__ import annotations

import time

from sqlalchemy import Column, Integer, String, Float, ForeignKey, Text, UniqueConstraint
from sqlalchemy.orm import relationship

from database.db import Base


class User(Base):
    __tablename__ = "users"

    id = Column(Integer, primary_key=True, index=True)
    email = Column(String(255), unique=True, index=True, nullable=False)
    name = Column(String(255), nullable=False, default="")
    password_hash = Column(String(255), nullable=False)
    role = Column(String(255), nullable=False, default="")
    institution = Column(String(255), nullable=False, default="")
    department = Column(String(255), nullable=False, default="")
    bio = Column(Text, nullable=False, default="")
    verified = Column(Integer, nullable=False, default=0)
    created_at = Column(Float, nullable=False, default=time.time)

    interests = relationship("ResearchInterest", back_populates="user", cascade="all, delete-orphan")
    saved_designs = relationship("SavedDesign", back_populates="user", cascade="all, delete-orphan")
    favorites = relationship("FavoriteGene", back_populates="user", cascade="all, delete-orphan")
    activity = relationship("RecentActivity", back_populates="user", cascade="all, delete-orphan")
    verification_tokens = relationship("VerificationToken", back_populates="user", cascade="all, delete-orphan")
    reports = relationship("Report", back_populates="user", cascade="all, delete-orphan")
    projects = relationship("Project", back_populates="user", cascade="all, delete-orphan")
    bug_reports = relationship("BugReport", back_populates="user", cascade="all, delete-orphan")


class ResearchInterest(Base):
    __tablename__ = "research_interests"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id"), nullable=False)
    topic = Column(String(255), nullable=False)
    description = Column(Text, nullable=False, default="")

    user = relationship("User", back_populates="interests")


class SavedDesign(Base):
    __tablename__ = "saved_designs"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id"), nullable=False)
    name = Column(String(255), nullable=False)
    gene_symbol = Column(String(100), nullable=False, default="")
    ensembl_id = Column(String(100), nullable=False, default="")
    disease = Column(String(255), nullable=False, default="")
    sequence = Column(Text, nullable=False, default="")
    notes = Column(Text, nullable=False, default="")
    created_at = Column(Float, nullable=False, default=time.time)

    user = relationship("User", back_populates="saved_designs")


class FavoriteGene(Base):
    __tablename__ = "favorite_genes"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id"), nullable=False)
    gene_symbol = Column(String(100), nullable=False)
    ensembl_id = Column(String(100), nullable=False, default="")
    note = Column(Text, nullable=False, default="")
    created_at = Column(Float, nullable=False, default=time.time)

    user = relationship("User", back_populates="favorites")


class RecentActivity(Base):
    __tablename__ = "recent_activity"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id"), nullable=False)
    action = Column(String(255), nullable=False)
    detail = Column(Text, nullable=False, default="")
    timestamp = Column(Float, nullable=False, default=time.time)

    user = relationship("User", back_populates="activity")


class Report(Base):
    __tablename__ = "reports"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id"), nullable=False)
    step = Column(String(50), nullable=False)
    title = Column(String(255), nullable=False)
    gene_symbol = Column(String(100), nullable=False, default="")
    disease = Column(String(255), nullable=False, default="")
    summary = Column(Text, nullable=False, default="")
    data = Column(Text, nullable=False, default="{}")
    created_at = Column(Float, nullable=False, default=time.time)

    user = relationship("User", back_populates="reports")


class VerificationToken(Base):
    __tablename__ = "verification_tokens"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id"), nullable=False)
    token = Column(String(255), unique=True, index=True, nullable=False)
    created_at = Column(Float, nullable=False, default=time.time)
    used = Column(Integer, nullable=False, default=0)

    user = relationship("User", back_populates="verification_tokens")


class Project(Base):
    __tablename__ = "projects"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id"), nullable=False)
    name = Column(String(255), nullable=False)
    description = Column(Text, nullable=False, default="")
    organism = Column(String(100), nullable=False, default="homo_sapiens")
    disease = Column(String(255), nullable=False, default="")
    gene_symbol = Column(String(100), nullable=False, default="")
    ensembl_id = Column(String(100), nullable=False, default="")
    therapeutic_goal = Column(String(255), nullable=False, default="")
    target_tissue = Column(String(255), nullable=False, default="")
    cell_line = Column(String(255), nullable=False, default="")
    notes = Column(Text, nullable=False, default="")
    status = Column(String(50), nullable=False, default="active")
    created_at = Column(Float, nullable=False, default=time.time)
    updated_at = Column(Float, nullable=False, default=time.time)

    user = relationship("User", back_populates="projects")


class BugReport(Base):
    __tablename__ = "bug_reports"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id"), nullable=False)
    area = Column(String(100), nullable=False)
    summary = Column(String(255), nullable=False)
    steps = Column(Text, nullable=False, default="")
    expected = Column(Text, nullable=False, default="")
    actual = Column(Text, nullable=False, default="")
    page_url = Column(String(500), nullable=False, default="")
    status = Column(String(50), nullable=False, default="open")
    created_at = Column(Float, nullable=False, default=time.time)

    user = relationship("User", back_populates="bug_reports")


class GeneFeatureBackup(Base):
    """Last-known-good gene structural feature analysis, keyed by gene.

    Used as a resilience backup so TG02 (Gene Function) keeps working for
    every gene when the Ensembl REST site is unreachable or a symbol cannot
    be resolved — the most recent successful analysis is replayed instead.
    """

    __tablename__ = "gene_feature_backups"

    id = Column(Integer, primary_key=True, index=True)
    organism = Column(String(100), nullable=False)
    gene_symbol = Column(String(100), nullable=False)
    ensembl_id = Column(String(100), nullable=False, default="")
    result = Column(Text, nullable=False, default="{}")
    created_at = Column(Float, nullable=False, default=time.time)
    updated_at = Column(Float, nullable=False, default=time.time)

    __table_args__ = (
        UniqueConstraint("organism", "gene_symbol", name="uq_gene_feature_backup"),
    )


class GeneLookupCache(Base):
    """Cached gene metadata keyed by organism + gene symbol.

    The full Ensembl expand=1 payload for a gene is expensive to build (e.g.
    DMD is ~300 KB and can take >15 s cold). Storing it means each gene is
    fetched from Ensembl once, then served instantly — and it doubles as a
    backup so a slow or unreachable Ensembl site never degrades the gene
    detail page for genes seen before.
    """

    __tablename__ = "gene_lookup_cache"

    id = Column(Integer, primary_key=True, index=True)
    organism = Column(String(100), nullable=False)
    gene_symbol = Column(String(100), nullable=False)
    payload = Column(Text, nullable=False, default="{}")
    created_at = Column(Float, nullable=False, default=time.time)
    updated_at = Column(Float, nullable=False, default=time.time)

    __table_args__ = (
        UniqueConstraint("organism", "gene_symbol", name="uq_gene_lookup_cache"),
    )


class RealDataCache(Base):
    """Last-known-good real data from an external source, keyed by namespace.

    Backs services/real_data_cache.py. The point is that an outage degrades to
    REAL DATA FETCHED EARLIER, or to an explicit "unavailable" — never to a
    synthesised substitute, which is indistinguishable from a measurement once
    it reaches the UI.

    `origin` separates a replayed live fetch from a hand-verified curated row;
    curated rows are never overwritten by a live fetch.
    """

    __tablename__ = "real_data_cache"

    id = Column(Integer, primary_key=True, index=True)
    namespace = Column(String(64), nullable=False, index=True)
    cache_key = Column(String(200), nullable=False, index=True)
    payload = Column(Text, nullable=False, default="{}")
    origin = Column(String(16), nullable=False, default="live")
    source = Column(String(200), nullable=False, default="")
    source_version = Column(String(64), nullable=False, default="")
    created_at = Column(Float, nullable=False, default=time.time)
    updated_at = Column(Float, nullable=False, default=time.time)

    __table_args__ = (
        UniqueConstraint("namespace", "cache_key", name="uq_real_data_cache"),
    )
